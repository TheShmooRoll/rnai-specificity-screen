import {decode, EXACT, labelTranscript, reverseComplement, transcriptCoordinates, WINDOW} from "./reference.js";
import {alignmentText, exactRuns, matches, scanReference, screenOptions, sharedKmerStarts} from "./matching.js";

const groupName = (transcript, gene) => transcript.gene === gene ? transcript.name : transcript.gene.name;
const groupId = (transcript, gene) => transcript.gene === gene ? transcript.id : transcript.gene.id;
const prefixCount = (array, start, end) => end > start ? array[end] - array[start] : 0;
const overlaps = (a, b) => a.start < b.end && b.start < a.end;

function targetSet(reference, gene, transcriptId) {
  if (!gene) return [];
  const isoforms = reference.isoforms(gene);
  if (!isoforms.length) throw new Error("This locus has no RNA model in the reference.");
  if (!transcriptId) return isoforms;
  const selected = isoforms.find(transcript => transcript.id === transcriptId);
  if (!selected) throw new Error("Select a transcript from this gene.");
  return [selected];
}

function selectInterval(hit, start, end) {
  return hit.start >= start && hit.end <= end;
}

function alignEvidence(reference, transcript, sequence, hit, origin) {
  const oriented = sequence.subarray(hit.start, hit.end);
  const insert = hit.reverse ? reverseComplement(oriented) : oriented;
  const subject = reference.sequence(transcript).subarray(hit.position, hit.position + hit.span);
  return {
    insertStart: hit.start - origin + 1, insertEnd: hit.end - origin,
    transcriptStart: hit.position + 1, transcriptEnd: hit.position + hit.span,
    orientation: hit.reverse ? "reverse complement" : "forward",
    edits: hit.edits, alignment: alignmentText(insert, subject, hit.gapped)
  };
}

function summarizeDirect(reference, gene, query, result, start, end) {
  const groups = new Map();
  for (const [type, size] of [["exact", EXACT], ["near", WINDOW]]) {
    for (const [target, strands] of result[type]) {
      const transcript = reference.transcripts[target];
      const id = groupId(transcript, gene);
      for (const hit of matches(strands, size, query.length)) {
        if (!selectInterval(hit, start, end)) continue;
        let group = groups.get(id);
        if (!group) {
          group = {id, name: groupName(transcript, gene), sibling: transcript.gene === gene,
            exactStarts: new Set(), nearStarts: new Set(), transcripts: new Map()};
          groups.set(id, group);
        }
        group[type + "Starts"].add(hit.start);
        let evidence = group.transcripts.get(target);
        if (!evidence) {
          evidence = {...labelTranscript(transcript), exactStarts: new Set(), nearStarts: new Set(),
            exactExample: null, nearExample: null};
          group.transcripts.set(target, evidence);
        }
        evidence[type + "Starts"].add(hit.start);
        const example = type + "Example";
        const previous = evidence[example];
        if (!previous || hit.edits < previous.edits) evidence[example] = hit;
      }
    }
  }
  const output = Array.from(groups.values(), group => ({
    id: group.id, name: group.name, sibling: group.sibling,
    exact: group.exactStarts.size, near: group.nearStarts.size,
    transcripts: Array.from(group.transcripts, ([target, value]) => ({
      id: value.id, name: value.name, type: value.type, length: value.length,
      exact: value.exactStarts.size, near: value.nearStarts.size,
      exactExample: value.exactExample && alignEvidence(reference, reference.transcripts[target], query, value.exactExample, start),
      nearExample: value.nearExample && alignEvidence(reference, reference.transcripts[target], query, value.nearExample, start)
    }))
  }));
  return output.sort((a, b) => b.exact - a.exact || b.near - a.near || a.name.localeCompare(b.name));
}

function matchingRuns(insert, target) {
  return [
    ...exactRuns(insert, target),
    ...exactRuns(reverseComplement(insert), target).map(run =>
      ({start: insert.length - run.end, end: insert.length - run.start, position: run.position}))
  ];
}

function mergedIntervals(intervals) {
  const merged = [];
  intervals.sort((a, b) => a.start - b.start || a.end - b.end);
  for (const interval of intervals) {
    const previous = merged.at(-1);
    if (previous && interval.start <= previous.end) {
      previous.end = Math.max(previous.end, interval.end);
    } else {
      merged.push({...interval});
    }
  }
  return merged;
}

function targetCoverage(reference, expected, insert, genomic) {
  return expected.map(transcript => {
    const target = reference.sequence(transcript);
    const runs = matchingRuns(insert, target);
    const insertRegions = mergedIntervals(runs.map(({start, end}) => ({start, end})));
    const starts = mergedIntervals(runs.map(({start, end}) => ({start, end: end - EXACT + 1})));
    const blocks = mergedIntervals(runs.map(run => ({start: run.position, end: run.position + run.end - run.start})));
    let exonicBases = null;
    if (genomic) {
      exonicBases = 0;
      for (let i = 0; i < transcript.exons.length; i += 2) {
        exonicBases += Math.max(0, Math.min(genomic.end, transcript.exons[i + 1]) -
          Math.max(genomic.start, transcript.exons[i]));
      }
    }
    return {...labelTranscript(transcript), exonicBases, blocks,
      matchedBases: insertRegions.reduce((sum, region) => sum + region.end - region.start, 0),
      exactStarts: starts.reduce((sum, region) => sum + region.end - region.start, 0),
      longestRun: runs.reduce((best, run) => Math.max(best, run.end - run.start), 0)
    };
  });
}

function publicCoverage(coverage) {
  return coverage.map(({blocks, ...model}) => ({...model,
    transcriptBlocks: blocks.map(block => ({start: block.start + 1, end: block.end}))}));
}

function summarizeUpstream(reference, gene, expected, coverage, screens) {
  const upstream = [];
  for (let i = 0; i < expected.length; i++) {
    const blocks = coverage[i].blocks;
    if (!blocks.length) continue;
    const transcript = expected[i];
    const sequence = reference.sequence(transcript);
    const nearest = new Map();
    for (const [type, size] of [["exact", EXACT], ["near", WINDOW]]) {
      for (const [target, strands] of screens[i][type]) {
        const subject = reference.transcripts[target];
        const id = groupId(subject, gene);
        for (const hit of matches(strands, size, sequence.length)) {
          const next = blocks.find(block => block.start >= hit.end);
          if (!next || blocks.some(block => overlaps(hit, block))) continue;
          const distance = next.start - hit.end;
          const previous = nearest.get(id);
          if (!previous || distance < previous.distance ||
            (distance === previous.distance && type === "exact" && previous.evidence !== "exact 21-base match")) {
            nearest.set(id, {id, name: groupName(subject, gene), transcript: subject.name,
              from: transcript.name, fromId: transcript.id, distance, triggerStart: next.start + 1,
              evidence: type === "exact" ? "exact 21-base match" : "40-base similarity",
              targetStart: hit.start + 1, targetEnd: hit.end});
          }
        }
      }
    }
    upstream.push(...nearest.values());
  }
  return upstream.sort((a, b) => a.distance - b.distance || a.name.localeCompare(b.name));
}

function genomicInterval(gene, start, end) {
  return gene.strand === "-" ? {start: gene.end - end, end: gene.end - start} :
    {start: gene.start + start, end: gene.start + end};
}

function prepareEvents(reference, gene, query, result) {
  const byGroup = new Map();
  const events = {exact: [], near: [], count: 0, byGroup, names: []};
  for (const [type, size] of [["exact", EXACT], ["near", WINDOW]]) {
    const unique = new Set();
    for (const [index, strands] of result[type]) {
      const transcript = reference.transcripts[index];
      const id = groupId(transcript, gene);
      if (!byGroup.has(id)) {
        byGroup.set(id, events.count++);
        events.names.push(groupName(transcript, gene));
      }
      const group = byGroup.get(id);
      for (const hit of matches(strands, size, query.length)) {
        const key = group * query.length + hit.start;
        if (unique.has(key)) continue;
        unique.add(key);
        (events[type][hit.start] ??= []).push(group);
      }
    }
  }
  return events;
}

function prepareSource(reference, gene, context, result, expected) {
  const query = context.sequence;
  const bad = new Uint32Array(query.length + 1);
  const coding = new Uint32Array(query.length + 1);
  const tx = context.transcript;
  for (let i = 0; i < query.length; i++) bad[i + 1] = bad[i] + (query[i] >= 4 ? 1 : 0);
  if (tx) {
    for (let i = 0; i < tx.coding.length; i += 2) {
      coding.fill(1, tx.coding[i] + 1, tx.coding[i + 1] + 1);
    }
    for (let i = 1; i < coding.length; i++) coding[i] += coding[i - 1];
  }
  return {...context, bad, coding, coordinates: tx && transcriptCoordinates(tx),
    coverage: expected.map(target => sharedKmerStarts(query, reference.sequence(target))),
    events: prepareEvents(reference, gene, query, result), result};
}

function visitCandidates(sources, min, max, visit, accept = null) {
  for (const source of sources) {
    const {events, sequence} = source;
    const counts = new Int32Array(events.count);
    for (let length = min; length <= Math.min(max, sequence.length); length++) {
      counts.fill(0);
      let groups = 0;
      let windows = 0;
      const exactCount = length - EXACT + 1;
      const nearCount = Math.max(0, length - WINDOW + 1);
      const update = (members, change) => {
        if (!members) return;
        for (const group of members) {
          if (counts[group] === 0) groups++;
          counts[group] += change;
          if (counts[group] === 0) groups--;
          windows += change;
        }
      };
      for (let i = 0; i < exactCount; i++) update(events.exact[i], 1);
      for (let i = 0; i < nearCount; i++) update(events.near[i], 1);
      for (let start = 0, end = length; end <= sequence.length; start++, end++) {
        if (source.bad[end] === source.bad[start] && (!accept || accept(source, start, end, counts, groups))) {
          let covered = 0;
          let shared = Infinity;
          for (const prefix of source.coverage) {
            const count = prefixCount(prefix, start, start + exactCount);
            if (count) covered++;
            shared = Math.min(shared, count);
          }
          if (covered) visit(source, start, end, covered, groups, windows, shared,
            source.coding[end] - source.coding[start]);
        }
        update(events.exact[start], -1);
        update(events.exact[start + exactCount], 1);
        if (nearCount) {
          update(events.near[start], -1);
          update(events.near[start + nearCount], 1);
        }
      }
    }
  }
}

function better(candidate, previous) {
  if (!previous) return true;
  return (previous.covered - candidate.covered || candidate.groups - previous.groups ||
    candidate.windows - previous.windows || previous.length - candidate.length ||
    previous.shared - candidate.shared || previous.coding - candidate.coding ||
    candidate.start - previous.start) < 0;
}

function coordinatesFor(source, gene, start, end) {
  if (!source.transcript) {
    const interval = genomicInterval(gene, start, end);
    return [{start: interval.start, end: interval.end}];
  }
  if (!source.coordinates) return [];
  const coordinates = source.coordinates;
  const step = gene.strand === "-" ? -1 : 1;
  const ranges = [];
  let first = start;
  for (let i = start + 1; i <= end; i++) {
    if (i < end && coordinates[i] === coordinates[i - 1] + step) continue;
    ranges.push({start: Math.min(coordinates[first], coordinates[i - 1]),
      end: Math.max(coordinates[first], coordinates[i - 1]) + 1});
    first = i;
  }
  return ranges;
}

function findIndependent(sources, gene, primary, direct, min, max, expected) {
  if (primary.covered !== expected) return null;
  const regions = coordinatesFor(primary.source, gene, primary.start, primary.end);
  if (!regions.length) return null;
  const original = primary.source.sequence.subarray(primary.start, primary.end);
  const tests = new Map(sources.map(source => {
    const overlap = new Uint32Array(source.sequence.length + 1);
    for (let i = 0; i < source.sequence.length; i++) {
      const position = source.transcript ? source.coordinates?.[i] :
        gene.strand === "-" ? gene.end - 1 - i : gene.start + i;
      overlap[i + 1] = overlap[i] + Number(regions.some(region => position >= region.start && position < region.end));
    }
    return [source, {overlap, shared: sharedKmerStarts(source.sequence, original),
      blocked: direct.map(group => source.events.byGroup.get(group.id)).filter(id => id !== undefined)}];
  }));
  let second = null;
  visitCandidates(sources, min, max, (source, start, end, covered, groups, windows, shared, coding) => {
    if (covered !== expected) return;
    const candidate = {source, start, end, length: end - start, covered, groups, windows, shared, coding};
    if (better(candidate, second)) second = candidate;
  }, (source, start, end, counts) => {
    if (source.transcript && !source.coordinates) return false;
    const {overlap, shared, blocked} = tests.get(source);
    return overlap[start] === overlap[end] && prefixCount(shared, start, end - EXACT + 1) === 0 &&
      blocked.every(index => counts[index] === 0);
  });
  return second;
}

function frontierEntry(candidate) {
  const {source, start, end, length, covered, groups} = candidate;
  const names = new Set();
  for (const [type, size] of [["exact", EXACT], ["near", WINDOW]]) {
    for (let i = start; i + size <= end; i++) {
      for (const group of source.events[type][i] ?? []) names.add(group);
    }
  }
  return {length, covered, offTargets: groups,
    offTargetNames: Array.from(names, group => source.events.names[group]).sort(),
    start: start + 1, end, transcript: source.transcript?.name ?? null};
}

function geneLabel(gene) {
  return gene && {id: gene.id, name: gene.name, strand: gene.strand, locus: gene.locus};
}

function candidateResult(reference, gene, candidate, expected, upstreamScreens, precomputedDirect) {
  const {source, start, end} = candidate;
  const insert = source.sequence.subarray(start, end);
  const genomic = source.transcript ? null : genomicInterval(gene, start, end);
  const coverage = targetCoverage(reference, expected, insert, genomic);
  const locations = coordinatesFor(source, gene, start, end).map(interval =>
    ({start: interval.start + 1, end: interval.end}));
  return {start: start + 1, end, sourceLength: source.sequence.length, length: insert.length,
    transcript: source.transcript && labelTranscript(source.transcript),
    chromosome: reference.chromosomes[gene.chromosome][0],
    locations, sequence: decode(insert), coverage: publicCoverage(coverage),
    direct: precomputedDirect ?? summarizeDirect(reference, gene, source.sequence, source.result, start, end),
    upstream: upstreamScreens ? summarizeUpstream(reference, gene, expected, coverage, upstreamScreens) : [],
    missing: coverage.filter(model => model.exactStarts === 0).map(model => model.name)};
}

export function checkInsert(reference, insert, options, onProgress) {
  const {gene, transcriptId, upstream = true} = options;
  const expected = targetSet(reference, gene, transcriptId);
  const exclude = new Set(expected.map(transcript => transcript.index));
  const settings = screenOptions(options.differences === undefined ? 2 : options.differences, options.indels);
  const context = [{sequence: insert, skipExcluded: true}];
  if (upstream) {
    for (const transcript of expected) context.push({sequence: reference.sequence(transcript), skipExcluded: true});
  }
  const [screen, ...upstreamScreens] = scanReference(reference, context, settings, exclude, onProgress);
  const coverage = targetCoverage(reference, expected, insert, null);
  return {mode: "check", gene: geneLabel(gene), settings, length: insert.length, coverage: publicCoverage(coverage),
    direct: summarizeDirect(reference, gene, insert, screen, 0, insert.length),
    upstream: upstream ? summarizeUpstream(reference, gene, expected, coverage, upstreamScreens) : []};
}

export async function designFragments(reference, options, onProgress = () => {}) {
  const {gene, transcriptId, min, max, template, upstream = true, findPair = true, findShorter = true} = options;
  if (!gene) throw new Error("Choose a target gene.");
  if (template !== "cdna" && template !== "genomic") throw new Error("Choose a cDNA or genomic DNA template.");
  if (!Number.isInteger(min) || !Number.isInteger(max) || min < EXACT || max < min || max > 5000) {
    throw new Error("Enter a valid fragment length or range between 21 and 5,000 bases.");
  }
  const expected = targetSet(reference, gene, transcriptId);
  const exclude = new Set(expected.map(transcript => transcript.index));
  const settings = screenOptions(options.differences === undefined ? 2 : options.differences, options.indels);
  const sourceContexts = template === "cdna" ? expected.filter(transcript => transcript.length >= EXACT)
    .map(transcript => ({sequence: reference.sequence(transcript), transcript, skipExcluded: true})) :
    [{sequence: await reference.genomic(gene), transcript: null, skipExcluded: true}];
  const info = {mode: "design", gene: geneLabel(gene), template, settings, expected: expected.map(labelTranscript)};
  const minFitsTemplate = sourceContexts.some(context => context.sequence.length >= min);
  const results = scanReference(reference, sourceContexts, settings, exclude, fraction => onProgress(fraction * .85));
  const sources = sourceContexts.map((context, i) => prepareSource(reference, gene, context, results[i], expected));
  const frontier = new Map();
  let primary = null;
  let partial = null;
  visitCandidates(sources, min, max, (source, start, end, covered, groups, windows, shared, coding) => {
    const length = end - start;
    const candidate = {source, start, end, length, covered, groups, windows, shared, coding};
    if (better(candidate, frontier.get(length))) frontier.set(length, candidate);
    if (better(candidate, primary)) primary = candidate;
    if (covered < expected.length && (!partial || groups < partial.groups ||
      (groups === partial.groups && better(candidate, partial)))) partial = candidate;
  });
  onProgress(.9);
  let shorter = null;
  if (findShorter && min > EXACT && (!primary || primary.groups || primary.covered < expected.length)) {
    visitCandidates(sources, EXACT, min - 1, (source, start, end, covered, groups, windows, shared, coding) => {
      if (covered !== expected.length) return;
      const candidate = {source, start, end, length: end - start, covered, groups, windows, shared, coding};
      if (better(candidate, shorter)) shorter = candidate;
    }, (_source, _start, _end, _counts, groups) => groups === 0);
  }
  const explanation = minFitsTemplate ? "No candidate within the requested range contains a contiguous 21-base target match." :
    "The requested minimum exceeds the length of every selected template.";
  const direct = primary && summarizeDirect(reference, gene, primary.source.sequence, primary.source.result, primary.start, primary.end);
  const independent = primary && findPair ? findIndependent(sources, gene, primary, direct, min, max, expected.length) : null;
  let upstreamScreens = null;
  if (upstream && (primary || shorter)) {
    upstreamScreens = expected.map(transcript => {
      const index = sourceContexts.findIndex(source => source.transcript === transcript);
      return index >= 0 ? results[index] : null;
    });
    const missing = expected.filter((_, i) => !upstreamScreens[i]);
    if (missing.length) {
      const scans = scanReference(reference, missing.map(transcript =>
        ({sequence: reference.sequence(transcript), skipExcluded: true})), settings, exclude,
      fraction => onProgress(.9 + fraction * .08));
      upstreamScreens = upstreamScreens.map(value => value ?? scans.shift());
    }
  }
  onProgress(1);
  return {...info, primary: primary && candidateResult(reference, gene, primary, expected, upstreamScreens, direct),
    independent: independent && candidateResult(reference, gene, independent, expected, upstreamScreens),
    alternative: partial && partial.groups < primary?.groups
      ? candidateResult(reference, gene, partial, expected, upstreamScreens) : null,
    shorter: shorter && candidateResult(reference, gene, shorter, expected, upstreamScreens),
    explanation: primary ? null : explanation,
    frontier: Array.from(frontier.values(), frontierEntry).sort((a, b) => a.length - b.length)};
}
