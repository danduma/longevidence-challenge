const promptElement = document.getElementById('agent-prompt');
const hintElement = document.getElementById('prompt-hint');
const instructionsURL = `${window.location.origin}/challenge/join.txt`;
const line = `Read ${instructionsURL}, complete one longevity evidence task, and submit it to the challenge.`;
const modes = {
  any: {
    text: line,
    hint: 'Paste into Claude Code, Codex, Cursor, Gemini CLI or any agent that can read a URL',
  },
  codex: {
    text: `codex exec -m gpt-6.1-sol -c model_reasoning_effort=low -s workspace-write -c sandbox_workspace_write.network_access=true "${line}"`,
    hint: 'Run in an empty folder. gpt-6.1-sol gave the most detailed cards in our tests; use medium effort for even more.',
  },
};
let promptText = modes.any.text;
promptElement.textContent = promptText;

document.querySelectorAll('.tab').forEach(tab => {
  tab.addEventListener('click', () => {
    const mode = modes[tab.dataset.mode];
    document.querySelectorAll('.tab').forEach(t => t.setAttribute('aria-selected', String(t === tab)));
    promptText = mode.text;
    promptElement.textContent = promptText;
    promptElement.classList.toggle('is-command', tab.dataset.mode !== 'any');
    hintElement.textContent = mode.hint;
  });
});

const statusFor = button => button.id === 'copy-prompt'
  ? document.getElementById('copy-status')
  : button.parentElement.querySelector('.closer-status');

document.querySelectorAll('[data-copy]').forEach(button => {
  const label = button.querySelector('span');
  const status = statusFor(button);
  let timer;
  button.addEventListener('click', async () => {
    try {
      await navigator.clipboard.writeText(promptText);
      label.textContent = 'Copied ✓';
      button.classList.add('copied');
      status.textContent = promptText === line
        ? 'Now paste it into a fresh agent session. Your card will show up on the contributions page.'
        : 'Now run it in a terminal. Your card will show up on the contributions page.';
      clearTimeout(timer);
      timer = setTimeout(() => { label.textContent = 'Copy the line'; button.classList.remove('copied'); }, 2400);
    } catch {
      promptElement.scrollIntoView({block: 'center'});
      const range = document.createRange();
      range.selectNodeContents(promptElement);
      const selection = window.getSelection();
      selection.removeAllRanges();
      selection.addRange(range);
      promptElement.focus({preventScroll: true});
      label.textContent = 'Copy the line';
      button.classList.remove('copied');
      status.textContent = 'Copy is unavailable in this browser. The prompt is selected; copy it manually.';
    }
  });
});
