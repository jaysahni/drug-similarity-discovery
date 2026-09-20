import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {test} from 'node:test';
const source = await readFile(new URL('../web/chats.js', import.meta.url), 'utf8');
const {createWorkspace, newChat, addFolder, updateChat, removeChat, removeFolder, findSavedResult, saveResult, removeSavedResult, readWorkspace, writeWorkspace} = await import(`data:text/javascript;base64,${Buffer.from(source).toString('base64')}`);
function memoryStorage() {
  const data = new Map();
  return {getItem: key => data.get(key) ?? null, setItem: (key, value) => data.set(key, value)};
}
test('new chats retain independent messages, drafts and result selections through reload', () => {
  const workspace = createWorkspace();
  const first = workspace.chats[0];
  Object.assign(first, {draft: 'unfinished question', messages: [{kind:'prompt', text:'Colorectal cancer'}], snapshotOpen: true, activeStep:'signature', viewState:{candidate:'42Q',view:'binding_site'}});
  const second = newChat(workspace);
  second.messages.push({kind:'prompt',text:'A different question'});
  assert.notEqual(first.id,second.id);
  assert.equal(second.draft,'');
  assert.equal(second.snapshotOpen,false);
  const storage = memoryStorage();
  assert.equal(writeWorkspace(storage,'test',workspace),true);
  const restored = readWorkspace(storage,'test');
  assert.equal(restored.writable,true);
  assert.deepEqual(restored.workspace,workspace);
  assert.equal(restored.workspace.chats[0].draft,'unfinished question');
  assert.equal(restored.workspace.chats[0].viewState.candidate,'42Q');
});
test('chats can be renamed, moved into folders and moved back without data loss', () => {
  const workspace = createWorkspace();
  const chat = workspace.chats[0];
  chat.messages.push({kind:'prompt',text:'Keep this message'});
  const folder = addFolder(workspace,' Oncology ');
  folder.collapsed = true;
  updateChat(workspace,chat.id,' My research ',folder.id);
  assert.equal(folder.name,'Oncology');
  assert.equal(folder.collapsed,false);
  assert.equal(chat.name,'My research');
  assert.equal(chat.folderId,folder.id);
  const nested = newChat(workspace,folder.id);
  removeFolder(workspace,folder.id);
  assert.equal(chat.folderId,null);
  assert.equal(nested.folderId,null);
  assert.equal(chat.messages[0].text,'Keep this message');
  assert.throws(()=>updateChat(workspace,chat.id,'Valid name','missing'));
  assert.throws(()=>addFolder(workspace,'   '));
});
test('deleting chats selects a survivor and always leaves a working tab', () => {
  const workspace = createWorkspace();
  const first = workspace.chats[0];
  const second = newChat(workspace);
  removeChat(workspace,second.id);
  assert.equal(workspace.activeId,first.id);
  removeChat(workspace,first.id);
  assert.equal(workspace.chats.length,1);
  assert.notEqual(workspace.activeId,first.id);
});
test('unreadable storage is preserved and denied writes are reported', () => {
  const storage = memoryStorage();
  storage.setItem('test','{broken');
  const restored=readWorkspace(storage,'test');
  assert.equal(restored.writable,false);
  assert.equal(restored.workspace.chats.length,1);
  assert.equal(storage.getItem('test'),'{broken');
  assert.equal(writeWorkspace({setItem(){throw Error('Quota exceeded');}},'test',restored.workspace),false);
  assert.equal(readWorkspace(null,'test').writable,false);
});
test('invalid saved relationships are rejected and run storage is isolated', () => {
  const storage=memoryStorage();
  const workspace=createWorkspace();
  writeWorkspace(storage,'run-a',workspace);
  assert.notEqual(readWorkspace(storage,'run-b').workspace.activeId,workspace.activeId);
  workspace.chats[0].folderId='missing';
  writeWorkspace(storage,'run-a',workspace);
  assert.equal(readWorkspace(storage,'run-a').writable,false);
});

test('saved results are independent bookmarks with deduplication and stable view state', () => {
  const workspace = createWorkspace();
  const chatId = workspace.activeId;
  const viewState = {candidate:'42Q', view:'binding_site'};
  const input = {runId:'3vhe-reference', candidateId:'42Q', name:'VEGFR2 binding example', subtitle:'42Q', chatId, viewState, activeStep:'signature'};
  const result = saveResult(workspace,input);
  viewState.view = 'candidate';
  assert.equal(result.viewState.view,'binding_site');
  const again = saveResult(workspace,{...input, viewState:{candidate:'42Q',view:'binding_site'}});
  assert.equal(again.id,result.id);
  assert.equal(workspace.savedResults.length,1);
  saveResult(workspace,{...input,candidateId:'other',subtitle:'Other molecule'});
  assert.equal(workspace.savedResults.length,2);
  removeChat(workspace,chatId);
  assert.equal(findSavedResult(workspace,input.runId,'42Q').id,result.id);
  const storage=memoryStorage();
  writeWorkspace(storage,'test',workspace);
  const restored=readWorkspace(storage,'test');
  assert.equal(restored.writable,true);
  assert.deepEqual(restored.workspace.savedResults,workspace.savedResults);
  removeSavedResult(restored.workspace,result.id);
  assert.equal(restored.workspace.savedResults.length,1);
  assert.equal(restored.workspace.chats.length,1);
});
test('existing chats migrate to an empty saved-results library without changes', () => {
  const workspace=createWorkspace();
  workspace.chats[0].draft='Keep my draft';
  const originalChat=structuredClone(workspace.chats[0]);
  delete workspace.savedResults;
  const storage=memoryStorage();
  writeWorkspace(storage,'test',workspace);
  const restored=readWorkspace(storage,'test');
  assert.equal(restored.writable,true);
  assert.deepEqual(restored.workspace.savedResults,[]);
  assert.deepEqual(restored.workspace.chats[0],originalChat);
});
test('invalid saved-result navigation is rejected without overwriting stored data', () => {
  const workspace=createWorkspace();
  saveResult(workspace,{runId:'run',candidateId:null,name:'Result',subtitle:'Target',chatId:workspace.activeId,viewState:{},activeStep:'invalid'});
  const storage=memoryStorage();
  writeWorkspace(storage,'test',workspace);
  const raw=storage.getItem('test');
  assert.equal(readWorkspace(storage,'test').writable,false);
  assert.equal(storage.getItem('test'),raw);
});
test('worked-example messages and open examples survive reload; invalid example fields are rejected', () => {
  const workspace = createWorkspace();
  const chat = workspace.chats[0];
  chat.messages.push({kind:'example', text:'Colorectal cancer', exampleId:'kdr'});
  chat.exampleId = 'kdr'; chat.exampleStep = 'search';
  saveResult(workspace, {runId:'example:kdr', candidateId:null, name:'Colorectal cancer', subtitle:'VEGFR2 · Candidate search', chatId: chat.id, viewState:{step:'search'}, activeStep:'search'});
  const storage = memoryStorage();
  assert.equal(writeWorkspace(storage, 'k', workspace), true);
  const restored = readWorkspace(storage, 'k');
  assert.equal(restored.writable, true);
  assert.deepEqual(restored.workspace, workspace);
  const patches = [
    w => {w.chats[0].exampleId = 'not valid!';},
    w => {w.chats[0].exampleStep = 'nowhere';},
    w => {w.chats[0].messages[0].exampleId = 42;},
    w => {w.chats[0].messages[0].kind = 'mystery';},
    w => {w.savedResults[0].activeStep = 'nowhere';},
  ];
  for (const patch of patches) {
    const broken = structuredClone(workspace); patch(broken);
    storage.setItem('broken', JSON.stringify(broken));
    assert.equal(readWorkspace(storage, 'broken').writable, false);
  }
  // Workspaces saved before worked examples existed still load.
  const legacy = structuredClone(workspace);
  delete legacy.chats[0].exampleId; delete legacy.chats[0].exampleStep;
  legacy.chats[0].messages = [{kind:'prompt', text:'Diabetes'}]; legacy.savedResults = [];
  storage.setItem('legacy', JSON.stringify(legacy));
  assert.equal(readWorkspace(storage, 'legacy').writable, true);
});
