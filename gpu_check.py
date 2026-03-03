#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import platform

def run(cmd):
    print(f"\n$ {' '.join(cmd)}")
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
        print(out.strip())
        return 0, out
    except subprocess.CalledProcessError as e:
        print(e.output.strip())
        return e.returncode, e.output
    except FileNotFoundError:
        print("command not found")
        return 127, ""

def main():
    print("=== Environment ===")
    print("python:", sys.executable)
    print("version:", sys.version.split()[0])
    print("platform:", platform.platform())
    print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))

    print("\n=== NVIDIA tools ===")
    print("nvidia-smi path:", shutil.which("nvidia-smi"))
    rc, _ = run(["nvidia-smi"])
    if rc != 0:
        print("\nNo working nvidia-smi. That usually means:")
        print("- you are on a CPU-only machine/runtime, OR")
        print("- your VS Code kernel is local (not Colab GPU), OR")
        print("- NVIDIA drivers aren't installed in that environment.")
    print("\n=== nvcc ===")
    print("nvcc path:", shutil.which("nvcc"))
    run(["nvcc", "--version"])

    print("\n=== PyTorch ===")
    try:
        import torch
        print("torch:", torch.__version__)
        print("torch.cuda.is_available():", torch.cuda.is_available())
        print("torch.cuda.device_count():", torch.cuda.device_count())
        if torch.cuda.is_available():
            print("GPU 0:", torch.cuda.get_device_name(0))
            x = torch.randn(4096, 4096, device="cuda")
            y = torch.randn(4096, 4096, device="cuda")
            z = x @ y  # matmul on GPU
            torch.cuda.synchronize()
            print("PyTorch GPU matmul: OK")
        else:
            print("PyTorch sees no CUDA GPU.")
    except Exception as e:
        print("PyTorch test failed:", repr(e))

    print("\n=== TensorFlow (optional) ===")
    try:
        import tensorflow as tf
        gpus = tf.config.list_physical_devices("GPU")
        print("tensorflow:", tf.__version__)
        print("TF GPUs:", gpus)
        if gpus:
            with tf.device("/GPU:0"):
                a = tf.random.normal([2048, 2048])
                b = tf.random.normal([2048, 2048])
                c = tf.matmul(a, b)
            print("TensorFlow GPU matmul: OK")
        else:
            print("TensorFlow sees no GPU (this is fine if you only use PyTorch/CUDA).")
    except Exception as e:
        print("TensorFlow test skipped/failed:", repr(e))

if __name__ == "__main__":
    main()