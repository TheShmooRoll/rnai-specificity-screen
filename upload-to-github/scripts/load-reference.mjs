import {createHash} from "node:crypto";
import {readFile} from "node:fs/promises";
import {resolve} from "node:path";
import {fileURLToPath} from "node:url";
import {gunzipSync} from "node:zlib";
import {makeReference} from "../src/reference.js";

export const bundledData = fileURLToPath(new URL("../src/data/", import.meta.url));

export async function loadReference(directory = bundledData) {
  const manifest = JSON.parse(await readFile(resolve(directory, "reference.json"), "utf8"));
  const pending = new Map();
  const read = filename => {
    if (!pending.has(filename)) {
      const result = readFile(resolve(directory, filename)).then(content => {
        const asset = manifest.assets[filename];
        if (!asset || asset.bytes !== content.length ||
          asset.sha256 !== createHash("sha256").update(content).digest("hex")) {
          throw new Error(`${filename} differs from the pinned reference manifest.`);
        }
        return gunzipSync(content);
      });
      pending.set(filename, result);
    }
    return pending.get(filename);
  };
  const [catalog, bases] = await Promise.all([read("catalog.json.gz"), read("transcripts.bin.gz")]);
  const reference = makeReference(JSON.parse(catalog.toString("utf8")), bases, read);
  return {reference, manifest};
}
