import assert from "node:assert/strict";
import {readFile} from "node:fs/promises";
import {scanReference, screenOptions} from "../src/matching.js";
import {EXACT, readInsert, reverseComplement} from "../src/reference.js";

const [release, negativeInsert, negativeTarget, positiveInsert, positiveTarget] = process.argv.slice(2);
if (!["WS235", "WS240"].includes(release) || !positiveTarget) {
  console.error("Usage: node scripts/check-clone-mapper.mjs WS235|WS240 negative-insert.fa negative-target.fa positive-insert.fa positive-target.fa");
  process.exit(1);
}

const checks = [
  {pair: "sjj_B0222.9 / F15E6.6", insert: negativeInsert, target: negativeTarget, positive: false},
  {pair: "sjj_Y50E8.g / ZK643.8a", insert: positiveInsert, target: positiveTarget, positive: true}
];
const results = [];
for (const check of checks) {
  const [insert, target] = await Promise.all([check.insert, check.target].map(async file =>
    readInsert(await readFile(file, "utf8"))));
  const transcript = {index: 0, length: target.length};
  const [screen] = scanReference({transcripts: [transcript], sequence: () => target},
    [{sequence: insert}], screenOptions(null));
  const observed = [];
  const expected = [];
  for (let strand = 0; strand < 2; strand++) {
    const query = strand ? reverseComplement(insert) : insert;
    const positions = [];
    for (let start = 0; start + EXACT <= query.length; start++) {
      const needle = query.subarray(start, start + EXACT);
      for (let offset = 0; offset + EXACT <= target.length; offset++) {
        if (needle.every((base, i) => base === target[offset + i])) {
          positions.push(start);
          break;
        }
      }
    }
    expected.push(positions);
    const matches = screen.exact.get(0)?.[strand];
    observed.push(matches ? [...matches.keys()].sort((a, b) => a - b) : []);
  }
  assert.deepEqual(observed, expected);
  const hits = observed[0].length + observed[1].length;
  assert.equal(hits > 0, check.positive, check.pair);
  // Use coordinates in the original insert for both orientations.
  results.push({pair: check.pair, forward: observed[0].map(start => start + 1),
    reverse: observed[1].map(start => insert.length - start - EXACT + 1).reverse()});
}
process.stdout.write(JSON.stringify({release, results}, null, 2) + "\n");
