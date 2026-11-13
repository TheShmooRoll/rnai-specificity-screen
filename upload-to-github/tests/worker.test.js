import assert from "node:assert/strict";
import test from "node:test";
import {readFile} from "node:fs/promises";

test("browser worker loads bundled data and completes the example insert workflow", async () => {
  const callbacks = new Map();
  const posted = [];
  const urls = [];
  const originalFetch = globalThis.fetch;
  globalThis.fetch = async url => {
    assert.equal(url.protocol, "file:");
    assert.match(url.pathname, /\/src\/data\//);
    urls.push(url.pathname.split("/").at(-1));
    return new Response(await readFile(url));
  };
  globalThis.self = {
    postMessage(message) {
      if (message.progress !== undefined) posted.push(message.progress);
      else callbacks.get(message.id)(message);
    }
  };
  try {
    await import("../src/worker.js");
    let id = 0;
    async function request(command, payload) {
      const current = ++id;
      const answer = new Promise(resolve => callbacks.set(current, resolve));
      await self.onmessage({data: {id: current, command, payload}});
      const {error, result} = await answer;
      assert.equal(error, undefined);
      return result;
    }
    const [gene] = await request("search", "kin-1");
    assert.equal(gene.name, "kin-1");
    assert.equal(gene.isoforms.length, 14);
    const example = await request("example", {gene: "ctl-2"});
    const result = await request("check", {gene: "ctl-2", insert: example.sequence, differences: 2, indels: true});
    assert.equal(result.reference.annotation, "WormBase WS298");
    assert.equal(result.settings.indels, true);
    assert.deepEqual(result.direct.slice(0, 2).map(hit => hit.name), ["ctl-1", "ctl-3"]);
    assert.deepEqual(urls.sort(), ["catalog.json.gz", "transcripts.bin.gz"]);
    assert.equal(posted[0], 0);
    assert.equal(posted.at(-1), 1);
  } finally {
    globalThis.fetch = originalFetch;
    delete globalThis.self;
  }
});
