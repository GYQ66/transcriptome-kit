# Agent 使用说明

本目录是一个 Agent Skill Suite（技能套件）。若你的客户端支持 Agent Skills
（如 WorkBuddy、Claude Code），把 `skills/` 下的各子目录复制到技能目录即可：

```bash
# WorkBuddy
cp -r skills/* ~/.workbuddy/skills/

# Claude Code
cp -r skills/* ~/.claude/skills/
```

复制后建议同时把根目录的 `SKILL.md` 也放进技能目录（可选，用于总控编排提示）：

```bash
cp SKILL.md ~/.workbuddy/skills/transcriptome-kit/
```

技能清单、依赖安装见 README.md 与 INSTALL.md。

## Agent 必读：环境检查与 R 安装询问

**首次在用户机器上使用本套件前，必须先运行：**

```bash
python <套件根目录>/check_env.py
```

- 检测到 R → 正常继续。
- **未检测到 R** → 脚本会打印 R 的下载地址并交互询问"是否自动下载"。
  **Agent 必须把这个询问原样转达给用户，等用户明确答复，不得代答、不得默认下载。**
  用户同意后运行 `python check_env.py --yes --simulate-missing` 所示的同一交互流程
  （或直接 `--print-download-url` 把直链给用户）；用户拒绝则把地址留给用户自行安装。
- 下载完成后 Windows 上还会**再问一次**"是否启动静默安装"，同样必须用户确认。
- 装好 R 后重跑 `python check_env.py --check-packages` 补查 R 包依赖。
