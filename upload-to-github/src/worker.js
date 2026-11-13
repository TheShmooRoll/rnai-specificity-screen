import {checkInsert, designFragments} from "./pipeline.js";
import {decode, findGene, labelTranscript, makeReference, readInsert} from "./reference.js";

async function loadAsset(name) {
  const response = await fetch(new URL("./data/" + name, import.meta.url));
  if (!response.ok) throw new Error(`Could not load ${name} (${response.status}).`);
  const bytes = new Uint8Array(await response.arrayBuffer());
  const compressed = bytes[0] === 0x1f && bytes[1] === 0x8b;
  const decoded = compressed ? new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))) :
    new Response(bytes);
  return name.endsWith(".json.gz") ? decoded.json() : new Uint8Array(await decoded.arrayBuffer());
}

let catalogPromise;
let referencePromise;
const chromosomes = new Map();
function getCatalog() {
  return catalogPromise ??= loadAsset("catalog.json.gz").catch(error => {
    catalogPromise = null;
    throw error;
  });
}
function getReference() {
  return referencePromise ??= Promise.all([getCatalog(), loadAsset("transcripts.bin.gz")]).then(([catalog, bases]) =>
    makeReference(catalog, bases, filename => {
      if (!chromosomes.has(filename)) chromosomes.set(filename, loadAsset(filename).catch(error => {
        chromosomes.delete(filename);
        throw error;
      }));
      return chromosomes.get(filename);
    })).catch(error => {
    referencePromise = null;
    throw error;
  });
}

function geneInfo(catalog, index) {
  const [id, name, locus, biotype, , , , , first, count] = catalog.genes[index];
  const isoforms = catalog.transcripts.slice(first, first + count).map(([id, name, type, length]) =>
    ({id, name, type, length}));
  return {id, name, locus, biotype, isoforms};
}

self.onmessage = async event => {
  const {id, command, payload} = event.data;
  const progress = fraction => self.postMessage({id, progress: fraction});
  try {
    let result;
    if (command === "search") {
      const catalog = await getCatalog();
      const term = payload.toLowerCase();
      const exact = [];
      const prefixes = [];
      for (let i = 0; i < catalog.genes.length; i++) {
        const [id, name, locus, , , , , , , count] = catalog.genes[i];
        if (!count) continue;
        const labels = [name, locus, id].map(value => value.toLowerCase());
        if (labels.includes(term)) exact.push(i);
        else if (prefixes.length < 8 && labels.some(value => value.startsWith(term))) prefixes.push(i);
      }
      result = [...exact, ...prefixes].slice(0, 8).map(index => geneInfo(catalog, index));
    } else {
      const reference = await getReference();
      const gene = payload.gene && findGene(reference.genes, payload.gene);
      if (payload.gene && !gene) throw new Error(`Could not find “${payload.gene}” in WS298.`);
      if (command === "example") {
        const transcript = reference.isoforms(gene)[0];
        result = {transcript: labelTranscript(transcript), sequence: decode(reference.sequence(transcript))};
      } else {
        const options = {...payload, gene};
        progress(0);
        result = command === "design" ? await designFragments(reference, options, progress) :
          checkInsert(reference, readInsert(payload.insert), options, progress);
        result.reference = {annotation: reference.metadata.annotation, assembly: reference.metadata.assembly,
          accession: reference.metadata.accession, annotationDate: reference.metadata.annotationDate};
      }
    }
    self.postMessage({id, result});
  } catch (error) {
    self.postMessage({id, error: error.message});
  }
};
