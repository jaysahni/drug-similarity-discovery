const id = () => crypto.randomUUID();
const cleanName = value => String(value).trim().slice(0, 80);
export const stepIds = ['literature', 'target', 'signature', 'candidates', 'disease', 'protein', 'search', 'check'];
const exampleIdValid = value => value === null || value === undefined || (typeof value === 'string' && /^[a-z0-9-]{1,40}$/.test(value));

export function newChat(workspace, folderId = null) {
  if (folderId !== null && !workspace.folders.some(f => f.id === folderId)) throw new Error('Folder not found');
  const chat = {id: id(), name: 'New chat', named: false, folderId, messages: [], draft: '', snapshotOpen: false, viewState: null, activeStep: 'candidates', exampleId: null, exampleStep: 'disease'};
  workspace.chats.push(chat);
  workspace.activeId = chat.id;
  return chat;
}
export function createWorkspace() {
  const workspace = {version: 1, activeId: null, folders: [], chats: [], savedResults: []};
  newChat(workspace);
  return workspace;
}
export function addFolder(workspace, name) {
  name = cleanName(name);
  if (!name) throw new Error('Enter a folder name');
  const folder = {id: id(), name, collapsed: false};
  workspace.folders.push(folder);
  return folder;
}
export function updateChat(workspace, chatId, name, folderId) {
  const chat = workspace.chats.find(c => c.id === chatId);
  if (!chat) throw new Error('Chat not found');
  name = cleanName(name);
  if (!name) throw new Error('Enter a chat name');
  if (folderId !== null && !workspace.folders.some(f => f.id === folderId)) throw new Error('Folder not found');
  Object.assign(chat, {name, named: true, folderId});
  const folder = workspace.folders.find(f => f.id === folderId);
  if (folder) folder.collapsed = false;
}
export function removeChat(workspace, chatId) {
  const index = workspace.chats.findIndex(c => c.id === chatId);
  if (index < 0) return;
  workspace.chats.splice(index, 1);
  if (!workspace.chats.length) newChat(workspace);
  else if (workspace.activeId === chatId) workspace.activeId = workspace.chats[Math.min(index, workspace.chats.length - 1)].id;
}
export function removeFolder(workspace, folderId) {
  workspace.folders = workspace.folders.filter(f => f.id !== folderId);
  workspace.chats.forEach(c => {if (c.folderId === folderId) c.folderId = null;});
}
export function findSavedResult(workspace, runId, candidateId) {
  return workspace.savedResults.find(result => result.runId === runId && result.candidateId === candidateId);
}
export function saveResult(workspace, {runId, candidateId, name, subtitle, chatId, viewState, activeStep}) {
  let result = findSavedResult(workspace, runId, candidateId);
  if (!result) {
    result = {id: id()};
    workspace.savedResults.push(result);
  }
  Object.assign(result, {runId, candidateId, name: cleanName(name), subtitle: String(subtitle).slice(0, 160), chatId,
    viewState: structuredClone(viewState), activeStep});
  return result;
}
export function removeSavedResult(workspace, resultId) {
  workspace.savedResults = workspace.savedResults.filter(result => result.id !== resultId);
}

export function readWorkspace(storage, key) {
  try {
    const raw = storage.getItem(key);
    if (!raw) return {workspace: createWorkspace(), writable: true};
    const data = JSON.parse(raw);
    if (data.version !== 1 || !Array.isArray(data.chats) || !data.chats.length || !Array.isArray(data.folders)) throw new Error('Invalid workspace');
    const ids = new Set();
    const validId = value => typeof value === 'string' && /^[a-zA-Z0-9-]{1,80}$/.test(value);
    const nameValid = value => typeof value === 'string' && value.trim().length > 0 && value.length <= 80;
    for (const f of data.folders) {
      if (!validId(f.id) || ids.has(f.id) || !nameValid(f.name) || typeof f.collapsed !== 'boolean') throw new Error('Invalid folder');
      ids.add(f.id);
    }
    const folders = new Set(ids);
    for (const c of data.chats) {
      if (!validId(c.id) || ids.has(c.id) || !nameValid(c.name) || typeof c.named !== 'boolean' || (c.folderId !== null && !folders.has(c.folderId)) || typeof c.draft !== 'string' || c.draft.length > 2000 || typeof c.snapshotOpen !== 'boolean' || !Array.isArray(c.messages)) throw new Error('Invalid chat');
      if (!stepIds.includes(c.activeStep)) throw new Error('Invalid stage');
      if (!exampleIdValid(c.exampleId) || (c.exampleStep !== undefined && c.exampleStep !== null && !stepIds.includes(c.exampleStep))) throw new Error('Invalid example view');
      if (c.viewState !== null && (typeof c.viewState !== 'object' || Array.isArray(c.viewState))) throw new Error('Invalid view');
      if (c.messages.some(m => !m || !['prompt','snapshot','example'].includes(m.kind) || typeof m.text !== 'string' || m.text.length > 2000)) throw new Error('Invalid message');
      if (c.messages.some(m => m.kind === 'example' && (typeof m.exampleId !== 'string' || !exampleIdValid(m.exampleId)))) throw new Error('Invalid example message');
      ids.add(c.id);
    }
    // Existing chat workspaces predate the saved-results library.
    if (data.savedResults === undefined) data.savedResults = [];
    if (!Array.isArray(data.savedResults)) throw new Error('Invalid saved results');
    const resultKeys = new Set();
    for (const result of data.savedResults) {
      const key = JSON.stringify([result.runId, result.candidateId]);
      if (!validId(result.id) || ids.has(result.id) || !nameValid(result.name)
          || typeof result.runId !== 'string' || !result.runId
          || (result.candidateId !== null && typeof result.candidateId !== 'string')
          || typeof result.subtitle !== 'string' || result.subtitle.length > 160
          || !validId(result.chatId) || resultKeys.has(key)
          || !stepIds.includes(result.activeStep)
          || !result.viewState || typeof result.viewState !== 'object' || Array.isArray(result.viewState)) throw new Error('Invalid saved result');
      ids.add(result.id); resultKeys.add(key);
    }
    if (!data.chats.some(c => c.id === data.activeId)) data.activeId = data.chats[0].id;
    return {workspace: data, writable: true};
  } catch {
    // Keep unreadable saved data intact instead of overwriting it with an empty chat.
    return {workspace: createWorkspace(), writable: false};
  }
}
export function writeWorkspace(storage, key, workspace) {
  try {storage.setItem(key, JSON.stringify(workspace)); return true;} catch {return false;}
}
