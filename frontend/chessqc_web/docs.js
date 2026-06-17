/* CHESS-QC in-site documentation viewer. Fetches the Markdown in /docs and renders it
 * with marked, themed by style.css. The docs are first-party trusted content, so the
 * rendered HTML is inserted directly. Works on a static host (GitHub Pages) and under the
 * launcher's local server: the docs live two levels up, at ../../docs/<NAME>.md. */
"use strict";

const DOCS = [
  { file: "USER_MANUAL",        title: "User Manual" },
  { file: "TECHNICAL_REFERENCE", title: "Technical Reference" },
  { file: "EQUATIONS",          title: "Equation Reference" },
  { file: "VALIDATION_REPORT",  title: "Validation Report" },
  { file: "OUTDATED_APPS",      title: "Outdated Applications" },
];
const GH_BASE = "https://github.com/Coastal-Hazards-System/CHESS-Quick-Compute/blob/main/docs/";

const $ = (id) => document.getElementById(id);
const nav = $("docNav"), body = $("docBody"), ghLink = $("ghLink");

// build sidebar
const links = {};
for (const d of DOCS) {
  const a = document.createElement("a");
  a.textContent = d.title;
  a.href = `?doc=${d.file}`;
  a.addEventListener("click", (e) => { e.preventDefault(); show(d.file); });
  nav.appendChild(a);
  links[d.file] = a;
}

function currentFile() {
  const q = new URLSearchParams(location.search).get("doc");
  return DOCS.some((d) => d.file === q) ? q : DOCS[0].file;
}

async function show(file) {
  const meta = DOCS.find((d) => d.file === file) || DOCS[0];
  for (const f in links) links[f].classList.toggle("active", f === meta.file);
  ghLink.href = GH_BASE + meta.file + ".md";
  history.replaceState(null, "", `?doc=${meta.file}`);
  document.title = `CHESS-QC — ${meta.title}`;
  body.innerHTML = "Loading…";
  try {
    const res = await fetch(`../../docs/${meta.file}.md`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const md = await res.text();
    marked.setOptions({ gfm: true, breaks: false });
    body.innerHTML = marked.parse(md);
    body.scrollTop = 0;
    window.scrollTo(0, 0);
  } catch (e) {
    body.innerHTML = `<p class="doc-error">Could not load this document (${e.message}). ` +
      `You can read it <a href="${GH_BASE + meta.file}.md" target="_blank" rel="noopener">on GitHub</a>.</p>`;
  }
}

show(currentFile());
CHESSQCUI.wireToggle("modeToggle", "modeLbl");
