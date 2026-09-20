# Pipeline diagram

`Pipeline.jsx` is the React component mounted in the saved-results conversation.
It reads the manifest and calls the existing result navigation when a user opens
a stage. It does not start analysis jobs or infer completion from missing data.

To edit and rebuild, run from `visualization/`:

```sh
npm ci
npm run build
```

The build writes `web/pipeline.js` and its license file. Both are included in
Python visualization exports, so viewing an export needs neither Node nor a CDN.
After rebuilding, create a fresh visualization export to include the changes.
