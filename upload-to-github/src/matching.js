import {decode, EXACT, reverseComplement, WINDOW} from "./reference.js";

const UINT32 = 2 ** 32;
const FILTER_MASK = (1 << 20) - 1;
const HIGH_MASK = (1 << (2 * (EXACT - 16))) - 1;

export function screenOptions(differences, indels = true) {
  if (differences === null) return {differences, indels: false, seed: 0, offsets: []};
  if (![0, 1, 2, 3].includes(differences)) throw new Error("Choose Off or an edit limit between 0 and 3.");
  const seed = [EXACT, 20, 13, 10][differences];
  const offsets = Array.from({length: differences + 1}, (_, i) => i * seed);
  return {differences, indels: indels && differences > 0, seed, offsets};
}

function keyFor(low, high, size) {
  if (size <= 16) return (low & ((2 ** (2 * size)) - 1)) >>> 0;
  return (high & ((1 << (2 * (size - 16))) - 1)) * UINT32 + low;
}

export function forEachKmer(sequence, size, visit) {
  let low = 0;
  let high = 0;
  let clean = 0;
  for (let end = 0; end < sequence.length; end++) {
    const base = sequence[end];
    clean = base < 4 ? clean + 1 : 0;
    high = ((high << 2) | (low >>> 30)) & HIGH_MASK;
    low = ((low << 2) | base) >>> 0;
    if (clean >= size) visit(keyFor(low, high, size), end - size + 1);
  }
}

export function sharedKmerStarts(query, other) {
  const keys = new Set();
  forEachKmer(other, EXACT, key => keys.add(key));
  forEachKmer(reverseComplement(other), EXACT, key => keys.add(key));
  const starts = new Uint32Array(query.length + 1);
  forEachKmer(query, EXACT, (key, start) => {
    if (keys.has(key)) starts[start + 1] = 1;
  });
  for (let i = 1; i < starts.length; i++) starts[i] += starts[i - 1];
  return starts;
}

function addSeeds(index, filter, query, size) {
  forEachKmer(query.sequence, size, (key, position) => {
    filter[key & FILTER_MASK] = 1;
    let positions = index.get(key);
    if (!positions) index.set(key, positions = []);
    positions.push({query, position});
  });
}

function countMismatches(query, start, reference, position, limit) {
  if (position < 0 || position + WINDOW > reference.length) return limit + 1;
  let edits = 0;
  for (let i = 0; i < WINDOW; i++) {
    if (query[start + i] >= 4 || query[start + i] !== reference[position + i]) {
      if (++edits > limit) break;
    }
  }
  return edits;
}

// Align an entire query flank to a reference prefix growing away from the seed.
// Only cells with |i-j| <= limit can be part of an alignment within the limit.
function alignFlank(query, qEdge, qLength, reference, rEdge, limit, direction, scratch) {
  const available = direction === 1 ? reference.length - rEdge : rEdge + 1;
  const columns = Math.min(available, qLength + limit);
  const infinity = limit + 1;
  if (columns < qLength - limit) return null;
  let previous = scratch[0];
  let current = scratch[1];
  for (let j = 0; j <= Math.min(columns, limit); j++) previous[j] = j;
  for (let i = 1; i <= qLength; i++) {
    const first = Math.max(0, i - limit);
    const last = Math.min(columns, i + limit);
    let rowMinimum = infinity;
    if (first === 0) current[0] = i;
    if (first === 0) rowMinimum = i;
    const qBase = query[qEdge + direction * (i - 1)];
    for (let j = Math.max(1, first); j <= last; j++) {
      const rBase = reference[rEdge + direction * (j - 1)];
      const substitution = previous[j - 1] + (qBase < 4 && qBase === rBase ? 0 : 1);
      const insertion = j > first ? current[j - 1] + 1 : infinity;
      const deletion = j < i + limit ? previous[j] + 1 : infinity;
      current[j] = Math.min(substitution, insertion, deletion);
      rowMinimum = Math.min(rowMinimum, current[j]);
    }
    if (rowMinimum > limit) return null;
    const swap = previous;
    previous = current;
    current = swap;
  }
  let best = infinity;
  let consumed = 0;
  for (let length = Math.max(0, qLength - limit); length <= columns; length++) {
    const score = previous[length];
    if (score < best || (score === best && Math.abs(length - qLength) < Math.abs(consumed - qLength))) {
      best = score;
      consumed = length;
    }
  }
  return best <= limit ? {edits: best, length: consumed} : null;
}

function verify(query, start, seedStart, reference, seedPosition, settings, maximum, scratch) {
  const before = seedStart - start;
  const after = WINDOW - before - settings.seed;
  const expectedStart = seedPosition - before;
  const mismatches = countMismatches(query, start, reference, expectedStart, maximum);
  let best = mismatches <= maximum
    ? {position: expectedStart, span: WINDOW, edits: mismatches}
    : null;
  const limit = best ? best.edits - 1 : maximum;
  if (!settings.indels || limit < 1) return best;
  const left = alignFlank(query, seedStart - 1, before, reference, seedPosition - 1, limit, -1, scratch);
  if (!left) return best;
  const right = alignFlank(query, seedStart + settings.seed, after, reference, seedPosition + settings.seed,
    limit - left.edits, 1, scratch);
  if (right) {
    best = {position: seedPosition - left.length, span: left.length + settings.seed + right.length,
      edits: left.edits + right.edits, gapped: true};
  }
  return best;
}

function record(exactOrNear, query, target, start, match) {
  let groups = exactOrNear.get(target);
  if (!groups) exactOrNear.set(target, groups = [new Map(), new Map()]);
  groups[query.reverse ? 1 : 0].set(start, match);
}

function findNear(seeds, sequence, referencePosition, target, settings, scratch, excluded) {
  for (const {query, position} of seeds) {
    if (excluded && query.skipExcluded) continue;
    for (const offset of settings.offsets) {
      const start = position - offset;
      if (start < 0 || start + WINDOW > query.sequence.length) continue;
      const previous = query.result.near.get(target)?.[query.reverse ? 1 : 0].get(start);
      if (previous?.edits === 0) continue;
      const limit = previous ? previous.edits - 1 : settings.differences;
      const hit = verify(query.sequence, start, position, sequence, referencePosition, settings, limit, scratch);
      if (hit) record(query.result.near, query, target, start, hit);
    }
  }
}

// Keep one best hit per query window, orientation and transcript. Every seed
// position is indexed, then only partition offsets anchor a 40-base window.
export function scanReference(reference, contexts, settings, excluded = new Set(), onProgress = () => {}) {
  const exactIndex = new Map();
  const nearIndex = settings.seed ? new Map() : null;
  const exactFilter = new Uint8Array(FILTER_MASK + 1);
  const nearFilter = nearIndex && new Uint8Array(FILTER_MASK + 1);
  const scratch = [new Uint8Array(WINDOW + 4), new Uint8Array(WINDOW + 4)];
  const results = contexts.map(context => {
    const result = {exact: new Map(), near: new Map()};
    for (const reverse of [false, true]) {
      const query = {...context, sequence: reverse ? reverseComplement(context.sequence) : context.sequence,
        reverse, result};
      if (context.exact !== false) addSeeds(exactIndex, exactFilter, query, EXACT);
      if (nearIndex && context.near !== false && query.sequence.length >= WINDOW) {
        addSeeds(nearIndex, nearFilter, query, settings.seed);
      }
    }
    return result;
  });
  const seed = settings.seed;
  const lowMask = seed <= 16 ? (2 ** (2 * seed)) - 1 : FILTER_MASK;
  for (const transcript of reference.transcripts) {
    const target = transcript.index;
    if (target % 1024 === 0) onProgress(target / reference.transcripts.length);
    if (transcript.length < EXACT) continue;
    const sequence = reference.sequence(transcript);
    const isExcluded = excluded.has(target);
    let clean = 0;
    let high = 0;
    let low = 0;
    for (let end = 0; end < sequence.length; end++) {
      const base = sequence[end];
      clean = base < 4 ? clean + 1 : 0;
      high = ((high << 2) | (low >>> 30)) & HIGH_MASK;
      low = ((low << 2) | base) >>> 0;
      if (clean >= EXACT && exactFilter[low & FILTER_MASK]) {
        const hits = exactIndex.get(high * UINT32 + low);
        if (hits) {
          const position = end - EXACT + 1;
          for (const {query, position: start} of hits) {
            if (isExcluded && query.skipExcluded) continue;
            const previous = query.result.exact.get(target)?.[query.reverse ? 1 : 0];
            if (!previous?.has(start)) {
              record(query.result.exact, query, target, start, {position, span: EXACT, edits: 0});
            }
          }
        }
      }
      if (nearIndex && clean >= seed && nearFilter[low & lowMask & FILTER_MASK]) {
        const hits = nearIndex.get(keyFor(low, high, seed));
        if (hits) findNear(hits, sequence, end - seed + 1, target, settings, scratch, isExcluded);
      }
    }
  }
  onProgress(1);
  return results;
}

export function exactRuns(source, target) {
  const index = new Map();
  forEachKmer(target, EXACT, (key, start) => {
    let positions = index.get(key);
    if (!positions) index.set(key, positions = []);
    positions.push(start);
  });
  const runs = [];
  forEachKmer(source, EXACT, (key, start) => {
    for (const position of index.get(key) ?? []) {
      const previous = source[start - 1];
      if (start && position && previous < 4 && previous === target[position - 1]) continue;
      let length = EXACT;
      while (start + length < source.length && position + length < target.length &&
        source[start + length] < 4 && source[start + length] === target[position + length]) length++;
      runs.push({start, end: start + length, position});
    }
  });
  return runs;
}

export function* matches(groups, length, total) {
  for (let strand = 0; strand < 2; strand++) {
    for (const [oriented, match] of groups[strand]) {
      const start = strand ? total - oriented - length : oriented;
      yield {...match, start, end: start + length, reverse: Boolean(strand)};
    }
  }
}

export function alignmentText(query, reference, gapped) {
  const a = decode(query);
  const b = decode(reference);
  if (!gapped) {
    const marks = Array.from(a, (base, index) => base !== "N" && base === b[index] ? "|" : "·");
    return [a, marks.join(""), b];
  }
  const table = Array.from({length: a.length + 1}, () => new Uint8Array(b.length + 1));
  for (let j = 0; j <= b.length; j++) table[0][j] = j;
  for (let i = 1; i <= a.length; i++) {
    table[i][0] = i;
    for (let j = 1; j <= b.length; j++) {
      table[i][j] = Math.min(table[i - 1][j - 1] + (a[i - 1] !== b[j - 1] || a[i - 1] === "N"),
        table[i - 1][j] + 1, table[i][j - 1] + 1);
    }
  }
  const top = [];
  const marks = [];
  const bottom = [];
  let i = a.length;
  let j = b.length;
  while (i || j) {
    const same = a[i - 1] === b[j - 1] && a[i - 1] !== "N";
    if (i && j && table[i][j] === table[i - 1][j - 1] + Number(!same)) {
      top.push(a[--i]);
      bottom.push(b[--j]);
      marks.push(same ? "|" : "·");
    } else if (j && table[i][j] === table[i][j - 1] + 1) {
      top.push("-");
      bottom.push(b[--j]);
      marks.push(" ");
    } else {
      top.push(a[--i]);
      bottom.push("-");
      marks.push(" ");
    }
  }
  return [top, marks, bottom].map(row => row.reverse().join(""));
}
