import {build} from 'esbuild';
import {readFile, writeFile} from 'node:fs/promises';
import {fileURLToPath} from 'node:url';

const root = fileURLToPath(new URL('../', import.meta.url));
await build({
  absWorkingDir: root,
  entryPoints: ['frontend/Pipeline.jsx'],
  outfile: 'web/pipeline.js',
  bundle: true,
  format: 'esm',
  target: ['es2020'],
  minify: true,
  define: {'process.env.NODE_ENV': '"production"'},
  legalComments: 'eof',
});
const licenses = await Promise.all(['react', 'react-dom', 'scheduler'].map(async name => {
  const license = await readFile(new URL(`../node_modules/${name}/LICENSE`, import.meta.url), 'utf8');
  return `${name}\n${license}`;
}));
await writeFile(new URL('../web/react-LICENSE.txt', import.meta.url), licenses.join('\n\n'));
console.log('Built the local React pipeline component.');
