import {readFile} from "node:fs/promises";
import {parseArgs} from "node:util";
import {checkInsert, designFragments} from "../src/pipeline.js";
import {findGene, readInsert} from "../src/reference.js";
import {loadReference} from "./load-reference.mjs";

const usage = `RNAi Specificity Screen · C. elegans WS298

  node scripts/screen.mjs check --insert insert.fa [--gene ctl-2]
  node scripts/screen.mjs design --gene ctl-2 --template cdna [--min 500 --max 1200]
  node scripts/screen.mjs design --gene kin-1 --template genomic --length 600

Common options
  --transcript ID             Select one RefSeq or WormBase transcript
  --differences off|0|1|2|3   Maximum edits in 40 insert bases (default 2)
  --ungapped                  Count only substitutions in the 40-base screen
  --no-upstream               Skip the upstream similarity search
  --data DIRECTORY            Use another copy of the prepared reference
  --no-pair                   Skip the search for a second fragment
  --help                      Print this help

The result is written to stdout as JSON. Pass --insert - for stdin.
`;

try {
  const {values, positionals} = parseArgs({
    allowPositionals: true,
    options: {
      insert: {type: "string"}, gene: {type: "string"}, template: {type: "string"},
      transcript: {type: "string"}, differences: {type: "string", default: "2"},
      min: {type: "string", default: "500"}, max: {type: "string", default: "1200"},
      length: {type: "string"}, data: {type: "string"},
      ungapped: {type: "boolean"}, "no-upstream": {type: "boolean"}, "no-pair": {type: "boolean"},
      help: {type: "boolean"}
    }
  });
  const [mode] = positionals;
  if (values.help || !mode) {
    process.stdout.write(usage);
  } else {
    if (!["check", "design"].includes(mode) || positionals.length !== 1) throw new Error("Choose check or design.");
    if (mode === "check" && !values.insert) throw new Error("Provide --insert FILE or --insert -.");
    if (mode === "design" && !values.gene) throw new Error("Provide a target with --gene.");
    if (!["off", "0", "1", "2", "3"].includes(values.differences)) {
      throw new Error("Choose a difference limit of off, 0, 1, 2 or 3.");
    }
    const {reference} = await loadReference(values.data);
    const gene = values.gene && findGene(reference.genes, values.gene);
    if (values.gene && !gene) throw new Error(`Could not find ${values.gene} in WS298.`);
    if (values.transcript && !gene) throw new Error("--transcript requires --gene.");
    const transcriptId = values.transcript && reference.isoforms(gene).find(model =>
      [model.id.toLowerCase(), model.name.toLowerCase()].includes(values.transcript.toLowerCase()))?.id;
    if (values.transcript && !transcriptId) throw new Error("This transcript does not belong to the selected locus.");
    const options = {gene, transcriptId, differences: values.differences === "off" ? null : Number(values.differences),
      indels: !values.ungapped, upstream: !values["no-upstream"], findPair: !values["no-pair"]};
    let result;
    if (mode === "design") {
      result = await designFragments(reference, {...options, template: values.template,
        min: Number(values.length ?? values.min), max: Number(values.length ?? values.max)});
    } else {
      const source = values.insert === "-" ? "/dev/stdin" : values.insert;
      result = checkInsert(reference, readInsert(await readFile(source, "utf8")), options);
    }
    const {assembly, annotation, accession, annotationDate} = reference.metadata;
    process.stdout.write(JSON.stringify({reference: {assembly, annotation, accession, annotationDate}, ...result}, null, 2) + "\n");
  }
} catch (error) {
  console.error(error.message);
  process.exitCode = 1;
}
