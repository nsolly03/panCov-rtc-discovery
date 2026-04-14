#!/bin/bash
# ================================================================
# Launch parallel CGCP monitoring session
# Left pane: your NIC5 account
# Right pane: Prof. Twizere NIC5 account
# Usage: bash cgcp_tmux_parallel.sh
# ================================================================

SESSION="cgcp-parallel"

# Kill existing session if any
tmux kill-session -t $SESSION 2>/dev/null

# Create new session
tmux new-session -d -s $SESSION -x 220 -y 50

# Left pane: your account
tmux send-keys -t $SESSION "ssh nic5" Enter
tmux send-keys -t $SESSION "watch -n 60 squeue -u onsekuye" Enter

# Right pane: Twizere's account
tmux split-window -h -t $SESSION
tmux send-keys -t $SESSION "ssh nic5-twizere" Enter
tmux send-keys -t $SESSION "watch -n 60 squeue -u twizere" Enter

# Attach
tmux attach-session -t $SESSION
