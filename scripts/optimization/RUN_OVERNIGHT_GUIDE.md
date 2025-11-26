# Running Overnight Script Without Sleep Interruption

This guide provides step-by-step instructions to run the overnight optimization script without your laptop sleeping affecting the execution.

## Method 1: Using `caffeinate` (Recommended for macOS)

The `caffeinate` command prevents your Mac from sleeping while the script runs.

### Steps:

1. **Open Terminal** and navigate to your project directory:
   ```bash
   cd /Users/danielmiles/Documents/School/Senior/Research/Material-Simulations-Research
   ```

2. **Run the script with caffeinate**:
   ```bash
   caffeinate -i ./scripts/optimization/run_overnight.sh 50
   ```
   
   Options:
   - `-i`: Prevents the system from idle sleeping
   - `-d`: Prevents the display from sleeping
   - `-m`: Prevents the disk from idle sleeping
   - `-s`: Prevents the system from sleeping (when on AC power)
   - `-u`: Prevents the system from sleeping (when on battery)
   
   **Recommended combination**:
   ```bash
   caffeinate -dis ./scripts/optimization/run_overnight.sh 50
   ```
   This prevents idle sleep, display sleep, and disk sleep.

3. **To run in background** (so you can close Terminal):
   ```bash
   caffeinate -dis ./scripts/optimization/run_overnight.sh 50 > overnight_run.log 2>&1 &
   ```

4. **Check progress**:
   ```bash
   tail -f overnight_run.log
   ```
   Or check the timestamped log:
   ```bash
   tail -f scripts/optimization/comparison_output/overnight_run_*/overnight_run.log
   ```

## Method 2: Using `screen` (Detachable Session)

This allows you to detach and reattach to the session, even if you close Terminal.

### Steps:

1. **Install screen** (if not already installed):
   ```bash
   # macOS - screen is usually pre-installed
   # If not, install via Homebrew: brew install screen
   ```

2. **Start a screen session**:
   ```bash
   screen -S optimization
   ```

3. **Run the script**:
   ```bash
   cd /Users/danielmiles/Documents/School/Senior/Research/Material-Simulations-Research
   caffeinate -dis ./scripts/optimization/run_overnight.sh 50
   ```

4. **Detach from screen** (keeps running in background):
   - Press `Ctrl+A` then `D`

5. **Reattach later** to check progress:
   ```bash
   screen -r optimization
   ```

6. **List all screen sessions**:
   ```bash
   screen -ls
   ```

## Method 3: Using `tmux` (Alternative to screen)

Similar to screen but more modern.

### Steps:

1. **Install tmux** (if not already installed):
   ```bash
   brew install tmux
   ```

2. **Start a tmux session**:
   ```bash
   tmux new -s optimization
   ```

3. **Run the script**:
   ```bash
   cd /Users/danielmiles/Documents/School/Senior/Research/Material-Simulations-Research
   caffeinate -dis ./scripts/optimization/run_overnight.sh 50
   ```

4. **Detach from tmux**:
   - Press `Ctrl+B` then `D`

5. **Reattach later**:
   ```bash
   tmux attach -t optimization
   ```

## Method 4: System Settings (Manual Prevention)

You can also manually prevent sleep in System Settings:

1. **Open System Settings** → **Battery** (or **Energy Saver** on older macOS)
2. **Set "Prevent automatic sleeping on power adapter when display is off"** to ON
3. **Keep your laptop plugged in**
4. **Set display sleep to "Never"** (or a very long time)

**Note**: This method requires you to manually change settings back after the run.

## Recommended Complete Workflow

Here's the recommended complete workflow:

```bash
# 1. Navigate to project directory
cd /Users/danielmiles/Documents/School/Senior/Research/Material-Simulations-Research

# 2. Start screen session with caffeinate
screen -S optimization
caffeinate -dis ./scripts/optimization/run_overnight.sh 50

# 3. Detach (Ctrl+A, then D)
# Script continues running even if you close Terminal

# 4. Later, check progress
screen -r optimization

# Or check logs directly
tail -f scripts/optimization/comparison_output/overnight_run_*/overnight_run.log
```

## Checking Progress Without Reattaching

You can check progress without reattaching to the session:

```bash
# Check if the process is running
ps aux | grep "run_overnight.sh"

# Check latest log file
ls -lt scripts/optimization/comparison_output/overnight_run_*/overnight_run.log | head -1 | xargs tail -f

# Or use the progress check script
./scripts/optimization/check_progress.sh
```

## Stopping the Script

If you need to stop the script:

1. **If attached to screen/tmux**: Press `Ctrl+C`
2. **If running in background**: Find the process and kill it:
   ```bash
   pkill -f "run_overnight.sh"
   # Or more specifically:
   ps aux | grep "run_overnight.sh" | grep -v grep | awk '{print $2}' | xargs kill
   ```

## Important Notes

1. **Keep laptop plugged in**: Even with caffeinate, it's best to keep your laptop plugged in
2. **Close the lid carefully**: On macOS, closing the lid may still cause sleep. Use `caffeinate -dis` to prevent this
3. **Check disk space**: Long runs can generate significant output. Ensure you have enough disk space
4. **Battery considerations**: If running on battery, the script will stop when the battery dies. Always plug in for overnight runs

## Example: Full Overnight Run

For a full overnight run with 100 evaluations:

```bash
# Start screen session
screen -S overnight_optimization

# Run with caffeinate (prevents sleep)
caffeinate -dis ./scripts/optimization/run_overnight.sh 100

# Detach: Ctrl+A, then D
# Script will continue running overnight

# Next morning, reattach to see results
screen -r overnight_optimization
```

## Troubleshooting

### Script stops when laptop sleeps
- Make sure you're using `caffeinate -dis` (or at least `caffeinate -i`)
- Keep laptop plugged in
- Check System Settings → Battery to ensure sleep prevention is enabled

### Can't find the screen session
```bash
screen -ls  # List all sessions
screen -r <session_name>  # Reattach to specific session
```

### Want to see output in real-time
```bash
# Reattach to screen/tmux session
screen -r optimization

# Or tail the log file
tail -f scripts/optimization/comparison_output/overnight_run_*/overnight_run.log
```

