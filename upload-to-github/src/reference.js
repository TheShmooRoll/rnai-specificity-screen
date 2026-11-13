const ALPHABET = "ACGTN";
const COMPLEMENT = [3, 2, 1, 0, 4];

export const EXACT = 21;
export const WINDOW = 40;

export function encode(sequence) {
  return Uint8Array.from(sequence, base => {
    const index = ALPHABET.indexOf(base === "U" ? "T" : base);
    return index === -1 ? 4 : index;
  });
}

export function decode(sequence) {
  return Array.from(sequence, base => ALPHABET[base]).join("");
}

export function reverseComplement(sequence) {
  return Uint8Array.from(sequence, (_, i) => COMPLEMENT[sequence[sequence.length - 1 - i]]);
}

export function readInsert(text) {
  const lines = text.trim().split(/\r?\n/);
  if (lines[0]?.startsWith(">")) lines.shift();
  const bases = lines.join("").replace(/\s/g, "").toUpperCase();
  if (!bases || /[^ACGTU]/.test(bases)) {
    throw new Error("Enter one DNA or RNA sequence using A, C, G and T or U. A single FASTA header is accepted.");
  }
  if (bases.length < EXACT || bases.length > 20000) {
    throw new Error("Enter an insert between 21 and 20,000 bases.");
  }
  return encode(bases);
}

export function makeReference(catalog, bases, loadChromosome) {
  if (catalog.metadata.schema !== 1 || bases.length !== catalog.metadata.transcriptBases) {
    throw new Error("The transcript file does not match the reference catalog.");
  }
  const genes = catalog.genes.map(([id, name, locus, biotype, chromosome, start, end, strand, first, count], index) =>
    ({index, id, name, locus, biotype, chromosome, start, end, strand, first, count}));
  const transcripts = catalog.transcripts.map(([id, name, type, length, exons, coding, gene, offset], index) =>
    ({index, id, name, type, length, exons, coding, gene: genes[gene], offset}));
  return {
    metadata: catalog.metadata, chromosomes: catalog.chromosomes, genes, transcripts, bases,
    sequence(transcript) {
      return bases.subarray(transcript.offset, transcript.offset + transcript.length);
    },
    isoforms(gene) {
      return transcripts.slice(gene.first, gene.first + gene.count);
    },
    async genomic(gene) {
      const [, , filename, length] = catalog.chromosomes[gene.chromosome];
      const chromosome = await loadChromosome(filename);
      if (chromosome.length !== length) throw new Error("The chromosome file does not match the reference catalog.");
      const slice = chromosome.subarray(gene.start, gene.end);
      return gene.strand === "-" ? reverseComplement(slice) : slice;
    }
  };
}

export function findGene(genes, term) {
  const query = term.trim().toLowerCase();
  const matches = genes.filter(gene => [gene.id, gene.name, gene.locus]
    .some(value => value && value.toLowerCase() === query));
  if (matches.length > 1) throw new Error("Several loci match. Select a specific WormBase gene ID.");
  return matches[0];
}

export function transcriptCoordinates(transcript) {
  const coordinates = new Int32Array(transcript.length);
  const direction = transcript.gene.strand === "-" ? -1 : 1;
  let position = 0;
  for (let e = 0; e < transcript.exons.length; e += 2) {
    const left = transcript.exons[e];
    const right = transcript.exons[e + 1];
    if (position + right - left > coordinates.length) return null;
    let base = direction === 1 ? left : right - 1;
    for (let i = left; i < right; i++, base += direction) coordinates[position++] = base;
  }
  return position === coordinates.length ? coordinates : null;
}

export function labelTranscript(transcript) {
  return {id: transcript.id, name: transcript.name, type: transcript.type, length: transcript.length};
}
