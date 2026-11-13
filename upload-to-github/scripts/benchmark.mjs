import {parseArgs} from "node:util";
import {checkInsert, designFragments} from "../src/pipeline.js";
import {encode, findGene} from "../src/reference.js";
import {loadReference} from "./load-reference.mjs";

const {values} = parseArgs({options: {
  count: {type: "string", default: "200"},
  seed: {type: "string", default: "20261113"},
  gene: {type: "string", multiple: true},
  data: {type: "string"}
}});
const count = Number(values.count);
const seed = Number(values.seed);
if (!Number.isSafeInteger(count) || count < 1 || !Number.isSafeInteger(seed)) {
  throw new Error("Pass a positive integer sample size and an integer seed.");
}

function randomGenerator(initial) {
  let state = initial >>> 0;
  return () => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    return state / 2 ** 32;
  };
}

function sample(genes, random, count) {
  const copy = [...genes];
  for (let i = 0; i < Math.min(count, copy.length); i++) {
    const swap = i + Math.floor(random() * (copy.length - i));
    [copy[i], copy[swap]] = [copy[swap], copy[i]];
  }
  return copy.slice(0, count);
}

function key(result) {
  const candidate = result.primary;
  return candidate ? [candidate.transcript.id, candidate.start, candidate.end].join("/") : null;
}

function hitKeys(result) {
  return new Set(result.direct.map(group => group.id));
}

const {reference} = await loadReference(values.data);
const eligible = reference.genes.filter(gene => gene.biotype === "protein_coding" &&
  reference.isoforms(gene).some(model => model.length >= 500));
const genes = values.gene?.map(name => {
  const gene = findGene(reference.genes, name);
  if (!gene) throw new Error(`Unknown gene ${name}`);
  return gene;
}) ?? sample(eligible, randomGenerator(seed), count);
let ungappedChanges = 0;
let indelChanges = 0;
let extraUngapped = 0;
let extraGapped = 0;
let totalMilliseconds = [0, 0, 0];
const output = [];

for (const gene of genes) {
  const settings = [
    {differences: null, indels: false},
    {differences: 2, indels: false},
    {differences: 2, indels: true}
  ];
  const results = [];
  for (let i = 0; i < settings.length; i++) {
    const start = performance.now();
    results.push(await designFragments(reference, {gene, template: "cdna", min: 500, max: 1200,
      upstream: false, findPair: false, findShorter: false, ...settings[i]}));
    totalMilliseconds[i] += performance.now() - start;
  }
  const [exact, ungapped, gapped] = results.map(key);
  const [exactInsert, ungappedInsert, gappedInsert] = results.map(result => result.primary?.sequence ?? null);
  ungappedChanges += Number(exactInsert !== ungappedInsert);
  indelChanges += Number(ungappedInsert !== gappedInsert);
  // Measure hit yield on the same insert even when the recommendation changes.
  const fixed = results[0].primary;
  let addedUngapped = 0;
  let addedGapped = 0;
  if (fixed) {
    const insert = encode(fixed.sequence);
    const [a, b, c] = settings.map(settings =>
      hitKeys(checkInsert(reference, insert, {gene, upstream: false, ...settings})));
    addedUngapped = Array.from(b).filter(id => !a.has(id)).length;
    addedGapped = Array.from(c).filter(id => !b.has(id)).length;
  }
  extraUngapped += addedUngapped;
  extraGapped += addedGapped;
  output.push({gene: gene.name, id: gene.id, exact, ungapped, gapped, addedUngapped, addedGapped});
  if (output.length % 10 === 0 || output.length === genes.length) {
    console.error(`Benchmarked ${output.length}/${genes.length} loci`);
  }
}

process.stdout.write(JSON.stringify({reference: reference.metadata.annotation, seed,
  samplingFrame: values.gene ? "Specified loci" : "Protein-coding loci with a mature RNA at least 500 bases long",
  eligible: eligible.length, count: genes.length, ungappedChanges, indelChanges,
  ungappedRate: ungappedChanges / genes.length, indelRate: indelChanges / genes.length,
  extraUngapped, extraGapped, meanDesignMilliseconds: totalMilliseconds.map(time => time / genes.length),
  loci: output
}, null, 2) + "\n");
