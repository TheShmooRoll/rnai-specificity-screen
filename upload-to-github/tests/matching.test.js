import test from "node:test";
import assert from "node:assert/strict";
import {matches, scanReference, screenOptions} from "../src/matching.js";
import {encode, reverseComplement} from "../src/reference.js";

function reference(sequences) {
  return {transcripts: sequences.map((sequence, index) => ({index, length: sequence.length, sequence})),
    sequence: transcript => transcript.sequence};
}

function random(initial) {
  let state = initial;
  return length => Uint8Array.from({length}, () => {
    state = (Math.imul(1664525, state) + 1013904223) >>> 0;
    return state >>> 30;
  });
}

function editDistance(a, b, indels) {
  if (!indels) {
    return a.length === b.length ? a.reduce((count, base, i) => count + (base > 3 || base !== b[i]), 0) : Infinity;
  }
  let previous = Array.from({length: b.length + 1}, (_, index) => index);
  for (let i = 0; i < a.length; i++) {
    const current = [i + 1];
    for (let j = 0; j < b.length; j++) {
      current.push(Math.min(previous[j + 1] + 1, current[j] + 1,
        previous[j] + (a[i] > 3 || a[i] !== b[j])));
    }
    previous = current;
  }
  return previous[b.length];
}

function actual(query, subject, limit, indels) {
  const [screen] = scanReference(reference([subject]), [{sequence: query}], screenOptions(limit, indels));
  const groups = screen.near.get(0);
  const coordinates = new Map();
  if (groups) {
    for (const hit of matches(groups, 40, query.length)) {
      const oriented = hit.reverse ? reverseComplement(query.subarray(hit.start, hit.end)) :
        query.subarray(hit.start, hit.end);
      const referenceSlice = subject.subarray(hit.position, hit.position + hit.span);
      assert.equal(editDistance(oriented, referenceSlice, indels), hit.edits);
      coordinates.set(`${hit.reverse}/${hit.start}`, hit.edits);
    }
  }
  return coordinates;
}

function exhaustive(query, subject, limit, indels) {
  const scores = new Map();
  for (const reverse of [false, true]) {
    for (let start = 0; start + 40 <= query.length; start++) {
      const slice = query.subarray(start, start + 40);
      const oriented = reverse ? reverseComplement(slice) : slice;
      let best = limit + 1;
      const drift = indels ? limit : 0;
      for (let pos = 0; pos + 40 - drift <= subject.length; pos++) {
        for (let span = 40 - drift; span <= Math.min(40 + drift, subject.length - pos); span++) {
          best = Math.min(best, editDistance(oriented, subject.subarray(pos, pos + span), indels));
        }
      }
      if (best <= limit) scores.set(`${reverse}/${start}`, best);
    }
  }
  return scores;
}

test("an internal insertion and a substitution require gapped verification", () => {
  const insert = random(12345)(40);
  const subject = Uint8Array.from([...insert.subarray(0, 19), 2, ...insert.subarray(19)]);
  subject[30] = (subject[30] + 1) % 4;
  const [exact] = scanReference(reference([subject]), [{sequence: insert}], screenOptions(null));
  assert.equal(exact.exact.size, 0);
  assert.equal(actual(insert, subject, 2, false).size, 0);
  assert.equal(actual(insert, subject, 1, true).size, 0);
  assert.equal(actual(insert, subject, 2, true).get("false/0"), 2);
  assert.equal(actual(reverseComplement(insert), subject, 2, true).get("true/0"), 2);
});

test("three spaced substitutions or indels are detected only at setting three", () => {
  const insert = random(9)(40);
  const substitutions = insert.slice();
  for (const position of [9, 19, 29]) substitutions[position] = (substitutions[position] + 1) % 4;
  const indels = Uint8Array.from(Array.from(insert).flatMap((base, i) =>
    i % 10 === 9 && i < 30 ? [3, base] : [base]));
  assert.equal(actual(insert, substitutions, 2, false).size, 0);
  assert.equal(actual(insert, substitutions, 3, false).get("false/0"), 3);
  assert.equal(actual(insert, indels, 2, true).size, 0);
  assert.equal(actual(insert, indels, 3, true).get("false/0"), 3);
});

test("seeded results agree with an exhaustive matcher across positions, orientations and edit settings", () => {
  const rng = random(32);
  const spaced = rng(73);
  const repetitive = encode("CAGACTCAGACTCAGACTCAGACTCAGACTCAGACTCAGACTCAGACT");
  const cases = [
    [spaced, Uint8Array.from(Array.from(spaced).flatMap((base, i) => i % 17 === 16 ? [2, base] : [base]))],
    [repetitive, Uint8Array.from([1, ...repetitive.slice(2, 18), 0, ...repetitive.slice(18, 40), 3])],
    [rng(42), rng(45)]
  ];
  for (let i = 0; i < 18; i++) {
    const query = rng(40 + i % 4);
    const oriented = i % 2 ? reverseComplement(query) : query;
    const subject = [...rng(i % 6), ...oriented, ...rng(i % 4)];
    const positions = rng(3);
    for (let edit = 0; edit < i % 4; edit++) {
      const location = Math.min(subject.length - 1, 4 + edit * 12 + positions[edit]);
      const base = (subject[location] + 1) % 4;
      if (edit === 0 && i % 3 === 1) subject.splice(location, 0, base);
      else if (edit === 0 && i % 3 === 2) subject.splice(location, 1);
      else subject[location] = base;
    }
    cases.push([query, Uint8Array.from(subject)]);
  }
  for (const [query, subject] of cases) {
    for (const limit of [0, 1, 2, 3]) {
      for (const indels of [false, true]) {
        assert.deepEqual(actual(query, subject, limit, indels), exhaustive(query, subject, limit, indels));
      }
    }
  }
});

test("exact matches honor sequence boundaries and search both orientations", () => {
  const query = encode("AGTCTGCTCCTACGGGCTATA");
  const sequences = [query.subarray(0, 11), query.subarray(11), reverseComplement(query),
    Uint8Array.from([4, ...query, 4]), Uint8Array.from([...query.subarray(0, 10), 4, ...query.subarray(11)])];
  const [screen] = scanReference(reference(sequences), [{sequence: query}], screenOptions(null));
  assert.deepEqual(Array.from(screen.exact.keys()), [2, 3]);
  assert.deepEqual(screen.exact.get(2)[1].get(0), {position: 0, span: 21, edits: 0});
  assert.deepEqual(screen.exact.get(3)[0].get(0), {position: 1, span: 21, edits: 0});
});
