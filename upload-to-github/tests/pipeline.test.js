import assert from "node:assert/strict";
import test from "node:test";
import {loadReference} from "../scripts/load-reference.mjs";
import {checkInsert, designFragments} from "../src/pipeline.js";
import {encode, findGene, makeReference, readInsert} from "../src/reference.js";

const bundled = loadReference();

test("WS298 ctl-2 transcript flags ctl-1 and ctl-3, and design identifies a shorter unflagged region", async () => {
  const {reference} = await bundled;
  const gene = findGene(reference.genes, "Y54G11A.5");
  assert.equal(gene.name, "ctl-2");
  const [transcript] = reference.isoforms(gene);
  const result = checkInsert(reference, reference.sequence(transcript), {gene, upstream: false});
  assert.deepEqual(result.direct.slice(0, 2).map(group => group.name), ["ctl-1", "ctl-3"]);
  assert.ok(result.direct.every(group => group.name !== "ctl-2"));
  assert.equal(result.coverage[0].matchedBases, transcript.length);

  const design = await designFragments(reference, {gene, template: "cdna", min: 500, max: 1200,
    upstream: true, findPair: true});
  assert.equal(design.primary.length, 571);
  assert.equal(design.primary.coverage[0].exactStarts, 551);
  assert.equal(design.primary.direct.length, 2);
  assert.deepEqual(design.frontier.find(row => row.length === 571).offTargetNames,
    design.primary.direct.map(group => group.name).sort());
  assert.equal(design.shorter.length, 403);
  assert.equal(design.shorter.direct.length, 0);
});

test("WS298 kin-1 genomic design counts introns and checks sibling isoforms", async () => {
  const {reference} = await bundled;
  const gene = findGene(reference.genes, "kin-1");
  const isoforms = reference.isoforms(gene);
  assert.equal(isoforms.length, 14);
  const design = await designFragments(reference, {gene, transcriptId: isoforms[0].id,
    template: "genomic", min: 600, max: 600, upstream: false});
  const candidate = design.primary;
  assert.equal(candidate.length, 600);
  assert.equal(candidate.coverage.length, 1);
  assert.ok(candidate.coverage[0].exonicBases < candidate.length);
  assert.ok(candidate.direct.some(group => group.sibling));
  assert.equal(candidate.locations[0].end - candidate.locations[0].start + 1, 600);
});

function smallReference() {
  const prefix = "GCTCGCCCGTGTGCACCAACTGACTTGTTCCT";
  const segment = "AGTCTGCTCCTACGGGCTATACATGCTCATATTCCTG";
  const a = encode(prefix + segment);
  const b = encode(prefix + "ACGTAC");
  const c = encode(segment.slice(0, 23));
  const rows = [
    ["t1", "t1", "mRNA", a.length, [0, a.length], [0, a.length], 0, 0],
    ["t2", "t2", "mRNA", b.length, [0, b.length], [0, b.length], 1, a.length],
    ["t3", "t3", "mRNA", c.length, [0, c.length], [0, c.length], 2, a.length + b.length]
  ];
  const metadata = {schema: 1, transcriptBases: a.length + b.length + c.length};
  const genes = rows.map((row, i) => [`g${i + 1}`, `g${i + 1}`, "", "protein_coding", 0,
    0, row[3], "+", i, 1]);
  const reference = makeReference({metadata, genes, chromosomes: [["I", "test", "test", a.length]],
    transcripts: rows}, Uint8Array.from([...a, ...b, ...c]), async () => a);
  return {reference, a, prefix, segment};
}

test("upstream distance is measured to the nearest downstream matching block", () => {
  const {reference, segment} = smallReference();
  const result = checkInsert(reference, encode(segment), {gene: reference.genes[0], differences: null});
  assert.equal(result.upstream.find(row => row.name === "g2")?.distance, 0);
  assert.equal(result.direct.find(row => row.name === "g2"), undefined);
  assert.equal(result.direct.find(row => row.name === "g3").exact, 3);
});

test("exact length can return a clearly labeled shorter fragment", async () => {
  const {reference} = smallReference();
  const result = await designFragments(reference, {gene: reference.genes[0], template: "cdna",
    min: 80, max: 80, differences: null, upstream: false});
  assert.equal(result.primary, null);
  assert.ok(result.shorter.length < 80);
  assert.equal(result.shorter.direct.length, 0);
});

test("FASTA accepts RNA and rejects a second record", () => {
  assert.equal(readInsert(">one\nAGUCUGCUCCUACGGGCUAUA\n").length, 21);
  assert.throws(() => readInsert(">one\nAGTC\n>two\nTGAC"), /one DNA or RNA sequence/);
});
