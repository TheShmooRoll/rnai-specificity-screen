import {readInsert} from "./reference.js";

const $ = id => document.getElementById(id);
const form = $("workspace");
const number = new Intl.NumberFormat("en");
const html = value => String(value ?? "").replace(/[&<>"']/g, char =>
  ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"})[char]);
const choice = name => form.querySelector(`input[name="${name}"]:checked`)?.value;

let worker;
let requestId = 0;
let pending = new Map();
let selectedGene;
let searchVersion = 0;
let highlighted = -1;
let matches = [];
let currentRun;
let lastResult;

function resetWorker(message) {
  worker?.terminate();
  worker = null;
  for (const item of pending.values()) item.reject(new Error(message));
  pending = new Map();
}

function request(command, payload, onProgress) {
  if (!worker) {
    worker = new Worker(new URL("./worker.js", import.meta.url), {type: "module"});
    worker.onmessage = ({data}) => {
      const item = pending.get(data.id);
      if (!item) return;
      if (data.progress !== undefined) {
        item.onProgress?.(data.progress);
        return;
      }
      pending.delete(data.id);
      if (data.error) item.reject(new Error(data.error));
      else item.resolve(data.result);
    };
    worker.onerror = () => resetWorker("The analysis worker stopped. Try again or reduce the search range.");
  }
  const id = ++requestId;
  return new Promise((resolve, reject) => {
    pending.set(id, {resolve, reject, onProgress});
    worker.postMessage({id, command, payload});
  });
}

function showError(message) {
  $("error").textContent = message;
  $("error").hidden = !message;
}

function renderForm() {
  const design = choice("workflow") === "design";
  const exact = choice("length-mode") === "exact";
  const difference = choice("differences");
  $("design-panel").hidden = !design;
  $("insert-panel").hidden = design;
  $("gene-optional").hidden = design;
  $("length-range").hidden = exact;
  $("length-exact").hidden = !exact;
  $("indels").disabled = difference === "off" || difference === "0";
  const limit = Number(difference);
  $("difference-hint").textContent = difference === "off" ? "The 40-base screen is off." :
    limit === 0 ? "Require a perfect 40-base match." :
    `At most ${limit} ${$("indels").checked ? "substitutions or indel bases" : "substitutions"} per window.`;
  if (limit === 3) $("difference-hint").textContent += " This setting can take longer.";
  $("run").innerHTML = currentRun ? "Analyzing…" :
    `${design ? "Find fragments" : "Check insert"} <span aria-hidden="true">→</span>`;
}

function closeSuggestions() {
  $("suggestions").hidden = true;
  $("gene").setAttribute("aria-expanded", "false");
  $("gene").removeAttribute("aria-activedescendant");
  highlighted = -1;
}

function resetGene() {
  selectedGene = null;
  $("isoforms").innerHTML = '<option value="">All annotated isoforms</option>';
  $("isoforms").disabled = true;
  $("isoform-hint").textContent = "Select a gene to see its annotated transcripts.";
  $("gene-meta").textContent = "Examples";
}

function selectGene(gene) {
  selectedGene = gene;
  searchVersion++;
  $("gene").value = gene.name;
  $("gene-meta").textContent = [gene.id, gene.locus].filter(Boolean).join(" · ");
  $("isoforms").innerHTML = '<option value="">All annotated isoforms</option>';
  for (const model of gene.isoforms) {
    $("isoforms").add(new Option(`${model.name} · ${number.format(model.length)} nt`, model.id));
  }
  $("isoforms").disabled = gene.isoforms.length < 2;
  const count = gene.isoforms.length;
  $("isoform-hint").textContent = count === 1 ? "One annotated transcript in this reference." :
    `${count} annotated transcripts. Selecting one screens the others as potential off-targets.`;
  closeSuggestions();
}

function highlight(index) {
  highlighted = (index + matches.length) % matches.length;
  for (let i = 0; i < matches.length; i++) {
    $("suggestion-" + i).setAttribute("aria-selected", String(i === highlighted));
  }
  $("gene").setAttribute("aria-activedescendant", "suggestion-" + highlighted);
}

async function searchGene(term, autoSelect = false) {
  const version = ++searchVersion;
  closeSuggestions();
  if (term.length < 2) return;
  try {
    const response = await request("search", term);
    if (version !== searchVersion) return;
    const exact = response.find(gene => [gene.id, gene.name, gene.locus]
      .some(value => value.toLowerCase() === term.toLowerCase()));
    if (autoSelect && exact) {
      selectGene(exact);
      return;
    }
    matches = response;
    $("suggestions").innerHTML = matches.map((gene, i) =>
      `<li id="suggestion-${i}" role="option" aria-selected="false" data-index="${i}">
        <strong>${html(gene.name)}</strong><span>${html(gene.locus || gene.id)}</span></li>`).join("");
    $("suggestions").hidden = !matches.length;
    $("gene").setAttribute("aria-expanded", String(Boolean(matches.length)));
    if (!matches.length) $("gene-meta").textContent = "No matching annotated gene";
  } catch (error) {
    if (version === searchVersion) showError(error.message);
  }
}

function updateInsertSize() {
  const bases = $("insert").value.trim().replace(/^>[^\n]*(?:\n|$)/, "").replace(/\s/g, "").length;
  $("insert-size").textContent = `${number.format(bases)} bases`;
}

$("gene").addEventListener("input", () => {
  resetGene();
  void searchGene($("gene").value.trim());
});
$("gene").addEventListener("keydown", event => {
  if (event.key === "Escape") closeSuggestions();
  if ($("suggestions").hidden) return;
  if (event.key === "ArrowDown" || event.key === "ArrowUp") {
    event.preventDefault();
    highlight(highlighted < 0 && event.key === "ArrowUp" ? matches.length - 1 :
      highlighted + (event.key === "ArrowDown" ? 1 : -1));
  }
  if (event.key === "Enter" && highlighted >= 0) {
    event.preventDefault();
    selectGene(matches[highlighted]);
  }
});
$("suggestions").addEventListener("mousedown", event => event.preventDefault());
$("suggestions").addEventListener("click", event => {
  const option = event.target.closest("li[data-index]");
  if (option) selectGene(matches[Number(option.dataset.index)]);
});
document.addEventListener("click", event => {
  if (!event.target.closest(".combobox")) closeSuggestions();
});
$("gene").addEventListener("blur", () => {
  searchVersion++;
  closeSuggestions();
});
for (const button of document.querySelectorAll("[data-gene]")) {
  button.addEventListener("click", () => {
    resetGene();
    $("gene").value = button.dataset.gene;
    void searchGene(button.dataset.gene, true);
  });
}

$("insert").addEventListener("input", updateInsertSize);
$("insert-example").addEventListener("click", async () => {
  const button = $("insert-example");
  button.disabled = true;
  try {
    const {transcript, sequence} = await request("example", {gene: "ctl-2"});
    $("insert").value = `>${transcript.name}\n${sequence.match(/.{1,70}/g).join("\n")}`;
    updateInsertSize();
    $("gene").value = "ctl-2";
    await searchGene("ctl-2", true);
  } catch (error) {
    showError(error.message);
  } finally {
    button.disabled = false;
  }
});

form.addEventListener("change", renderForm);

function setBusy(busy) {
  $("run").disabled = busy;
  $("cancel").hidden = !busy;
  $("progress-row").hidden = !busy;
  renderForm();
}

form.addEventListener("submit", async event => {
  event.preventDefault();
  if (currentRun) return;
  showError("");
  closeSuggestions();
  const workflow = choice("workflow");
  const gene = selectedGene?.id ?? $("gene").value.trim();
  const template = choice("template");
  if (workflow === "design" && !gene) {
    showError("Enter a target gene.");
    $("gene").focus();
    return;
  }
  if (workflow === "design" && !template) {
    showError("Choose cDNA or genomic DNA.");
    form.querySelector('[name="template"]').focus();
    return;
  }
  if (workflow === "check") {
    try {
      readInsert($("insert").value);
    } catch (error) {
      showError(error.message);
      $("insert").focus();
      return;
    }
  }
  const difference = choice("differences");
  const exact = choice("length-mode") === "exact";
  const payload = {gene, template, transcriptId: $("isoforms").value, insert: $("insert").value,
    min: Number($(exact ? "exact-length" : "min").value),
    max: Number($(exact ? "exact-length" : "max").value),
    differences: difference === "off" ? null : Number(difference),
    indels: $("indels").checked, upstream: $("upstream").checked};
  const run = Symbol();
  currentRun = run;
  $("progress-label").textContent = "Loading WS298 reference…";
  $("progress").removeAttribute("value");
  setBusy(true);
  try {
    const output = await request(workflow, payload, fraction => {
      $("progress-label").textContent = fraction >= .85 && workflow === "design" ?
        "Ranking candidate fragments…" : "Screening the transcriptome…";
      $("progress").value = Math.floor(fraction * 100);
    });
    if (currentRun !== run) return;
    lastResult = output;
    renderResults(output);
    $("results-title").focus();
  } catch (error) {
    if (currentRun === run) showError(error.message);
  } finally {
    if (currentRun === run) {
      currentRun = null;
      setBusy(false);
    }
  }
});

$("cancel").addEventListener("click", () => {
  currentRun = null;
  resetWorker("Analysis canceled.");
  setBusy(false);
  $("run").focus();
});

function save(name, content, type) {
  const url = URL.createObjectURL(new Blob([content], {type}));
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

$("download").addEventListener("click", () => {
  if (lastResult) save("rnai-specificity-results.json", JSON.stringify(lastResult, null, 2), "application/json");
});

function table(headers, rows) {
  return `<div class="table-wrap"><table><thead><tr>${headers.map(name => `<th scope="col">${name}</th>`).join("")}</tr></thead>
    <tbody>${rows.map(row => `<tr>${row.map(value => `<td>${value}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`;
}

function coverageTable(coverage) {
  if (!coverage.length) return "";
  const genomic = coverage.some(row => row.exonicBases !== null);
  const headers = ["Target transcript", ...(genomic ? ["Exonic bases"] : []),
    "Exact-matching insert bases", "Exact 21-base starts", "Longest exact run"];
  return table(headers, coverage.map(row => [
    `<b>${html(row.name)}</b><small>${html(row.id)}</small>`,
    ...(genomic ? [number.format(row.exonicBases)] : []), number.format(row.matchedBases),
    number.format(row.exactStarts), `${number.format(row.longestRun)} nt`
  ]));
}

function evidenceMarkup(example, label) {
  return example ? `<div class="evidence"><b>${label}</b><small>
      Insert ${example.insertStart}–${example.insertEnd} · Transcript ${example.transcriptStart}–${example.transcriptEnd}
      · ${html(example.orientation)} · ${example.edits} difference(s)
    </small><pre>${example.alignment.map(html).join("\n")}</pre></div>` : "";
}

function directMarkup(groups, settings, showReferenceMatches) {
  const title = showReferenceMatches ? "Reference matches" : "Direct off-targets";
  if (!groups.length) return `<h4>${title}</h4><p class="result-note">No matches found under the selected rules.</p>`;
  return `<h4>${title} <span class="count">${number.format(groups.length)}</span></h4>
    <div class="hit-list">${groups.map(group => `<details><summary>
      <span class="hit-name">${html(group.name)} ${group.sibling ? "<small>Sibling isoform</small>" : ""}</span>
      <span class="hit-counts">${number.format(group.exact)} exact${settings.differences === null ? "" : ` · ${number.format(group.near)} near`}</span>
      </summary><div class="hit-evidence">${group.transcripts.map(model => `<div>
        <b>${html(model.name)}</b> <small>${html(model.id)} · ${number.format(model.exact)} exact${settings.differences === null ? "" : ` · ${number.format(model.near)} near`}</small>
        ${evidenceMarkup(model.exactExample, "Exact 21-base example")}
        ${evidenceMarkup(model.nearExample, "Best 40-base example")}
      </div>`).join("")}</div></details>`).join("")}</div>`;
}

function upstreamMarkup(upstream) {
  if (!upstream.length) return "";
  return `<details class="result-details"><summary>Shared sequence upstream <span class="count">${upstream.length}</span></summary>
    <p>Distances are measured along the mature target transcript.</p>
    ${table(["Target transcript", "Shared with", "Distance", "Evidence"],
      upstream.map(row => [html(row.from), `${html(row.name)}<small>${html(row.transcript)}</small>`,
        `${row.distance} nt`, html(row.evidence)]))}
  </details>`;
}

function candidateMarkup(candidate, title, key, result) {
  const covered = candidate.coverage.length - candidate.missing.length;
  const source = candidate.transcript ? html(candidate.transcript.name) : `Chr ${html(candidate.chromosome)}`;
  const location = candidate.transcript ? `${number.format(candidate.start)}–${number.format(candidate.end)}` :
    candidate.locations.map(range => `${number.format(range.start)}–${number.format(range.end)}`).join(", ");
  const left = (candidate.start - 1) / candidate.sourceLength * 100;
  const width = candidate.length / candidate.sourceLength * 100;
  return `<article class="candidate">
    <div class="candidate-head"><div><h3>${title}</h3><p><b>${source}</b> · ${location} · ${number.format(candidate.length)} bp</p></div>
      <button type="button" class="secondary" data-fasta="${key}">Download FASTA</button></div>
    <div class="fragment-track" aria-hidden="true"><span style="left:${left}%;width:${width}%"></span></div>
    <div class="candidate-stats">
      <span><b>${covered}/${candidate.coverage.length}</b> isoforms covered</span>
      <span><b>${candidate.direct.length}</b> direct flags</span>
      <span><b>${number.format(Math.min(...candidate.coverage.map(row => row.matchedBases)))}</b>
        ${candidate.coverage.length === 1 ? "exact-matching bases" : "matching bases in the least-covered isoform"}</span>
      ${candidate.length < 500 ? '<span class="caution">Short insert · review efficacy</span>' : ""}
    </div>
    ${candidate.missing.length ? `<p class="notice">Missing exact target match in ${candidate.missing.map(html).join(", ")}.</p>` : ""}
    <details class="result-details"><summary>Target coverage and sequence</summary>
      ${coverageTable(candidate.coverage)}
      <pre class="insert-sequence">${html(candidate.sequence.match(/.{1,70}/g).join("\n"))}</pre>
    </details>
    ${directMarkup(candidate.direct, result.settings, false)}
    ${upstreamMarkup(candidate.upstream)}
  </article>`;
}

function frontierMarkup(result) {
  const values = result.frontier;
  if (values.length < 2) return "";
  const runs = [];
  for (const item of values) {
    const kind = item.covered < result.expected.length ? "partial" : item.offTargets ? "flagged" : "clear";
    const run = runs.at(-1);
    if (run?.kind === kind) run.count++;
    else runs.push({kind, count: 1});
  }
  const index = values.findIndex(item => item.length === result.primary.length);
  return `<details class="result-details frontier"><summary>Compare lengths</summary>
    <div class="length-track" aria-hidden="true">${runs.map(run => `<span class="${run.kind}" style="flex:${run.count}"></span>`).join("")}</div>
    <div class="track-labels"><span>${values[0].length} bp</span><span>${values.at(-1).length} bp</span></div>
    <div class="track-legend"><span class="clear">Full coverage, no flags</span><span class="flagged">Full coverage, direct flags</span><span class="partial">Partial coverage</span></div>
    <label class="frontier-label">Inspect a length <input type="range" id="frontier-length" min="0" max="${values.length - 1}" value="${index}"></label>
    <output id="frontier-output" for="frontier-length"></output>
  </details>`;
}

function renderResults(result) {
  $("download").hidden = false;
  $("results-subtitle").textContent = [
    result.gene?.name,
    result.mode === "design" ? (result.template === "cdna" ? "cDNA" : "Genomic DNA") :
      `${number.format(result.length)}-base insert`,
    result.reference.annotation
  ].filter(Boolean).join(" · ");
  const content = $("results-content");
  if (result.mode === "check") {
    content.innerHTML = `<div class="candidate">${result.coverage.length ? "<h4>Target coverage</h4>" : ""}
      ${coverageTable(result.coverage)}
      ${directMarkup(result.direct, result.settings, !result.gene)}
      ${upstreamMarkup(result.upstream)}
    </div>`;
    return;
  }
  const hasClear = result.primary && !result.primary.missing.length && !result.primary.direct.length;
  const shorter = result.shorter ? candidateMarkup(result.shorter, "Shorter fragment outside the requested range", "shorter", result) : "";
  content.innerHTML = result.primary ? `${hasClear ? "" :
      '<p class="notice">No candidate in this range covered the requested isoforms and passed the direct screen.</p>'}
    ${candidateMarkup(result.primary, hasClear ? "Recommended fragment" : "Best available fragment", "primary", result)}
    ${result.independent ? candidateMarkup(result.independent, "Independent fragment", "independent", result) :
      '<p class="result-note pair-note">No second fragment in this range met the coverage, coordinate, 21-mer and shared off-target checks.</p>'}
    ${result.alternative ? candidateMarkup(result.alternative, "Partial-coverage alternative", "alternative", result) : ""}
    ${shorter}${frontierMarkup(result)}` :
    `<p class="notice">${html(result.explanation)}</p>${shorter}`;
  for (const button of content.querySelectorAll("[data-fasta]")) {
    const candidate = result[button.dataset.fasta];
    button.addEventListener("click", () => {
      const name = `${result.gene.name}-${button.dataset.fasta}`;
      save(`${name}.fa`, `>${name} ${result.reference.annotation}\n${candidate.sequence.match(/.{1,70}/g).join("\n")}\n`, "text/plain");
    });
  }
  const slider = $("frontier-length");
  if (slider) {
    const update = () => {
      const item = result.frontier[Number(slider.value)];
      const description = `${item.length} bp · ${item.covered}/${result.expected.length} isoforms · ${item.offTargets} direct flags`;
      slider.setAttribute("aria-valuetext", description);
      $("frontier-output").textContent = item.offTargetNames.length ?
        `${description} · ${item.offTargetNames.join(", ")}` : description;
    };
    slider.addEventListener("input", update);
    update();
  }
}

renderForm();
