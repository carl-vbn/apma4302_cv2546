import subprocess
import re
import numpy as np
import os


def run_bvp(m, gamma, k):
    # This is a little cursed but I didnt find
    # a good way to export all the data in one execution
    cmd = [
        "mpiexec", "-np", "4", "./bvp",
        "-options_file", "options_file",
        "-bvp_m", str(m),
        "-bvp_gamma", str(gamma),
        "-bvp_k", str(k),
        "-bvp_c", "0.0",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout + result.stderr
    match = re.search(r"relative error.*=\s*([0-9.eE+-]+)", output)
    if match:
        return float(match.group(1))
    else:
        raise RuntimeError(f"Could not parse error from output:\n{output}")


def plot_convergence(k_values, h_values, all_errors, all_orders):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))

    for k, errors, orders in zip(k_values, all_errors, all_orders):
        avg_order = np.mean(orders)
        ax.loglog(h_values, errors, 'o-', label=f'k={k} (order {avg_order:.2f})')

    # Reference slope
    ax.loglog(h_values, h_values**2 * (all_errors[0][0] / h_values[0]**2) * 0.5,
              'k--', alpha=0.4, label='$O(h^2)$')

    ax.set_xlabel('h', fontsize=14)
    ax.set_ylabel('Relative error', fontsize=14)
    ax.set_title(r'Convergence of BVP error ($\gamma=0$)', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, which='both', alpha=0.3)

    plt.show()


if __name__ == "__main__":
    gamma = 0.0
    k_values = [1, 5, 10]
    m_values = [40, 80, 160, 320, 640, 1280]

    h_values = np.array([1.0 / (m - 1) for m in m_values])

    # Make sure executable exists
    if not os.path.isfile("./bvp"):
        raise FileNotFoundError("Executable 'bvp' not found. Please compile it first.")

    all_errors = []
    all_orders = []
    for k in k_values:
        print(f"Running BVP for k={k}...")
        errors = np.array([run_bvp(m, gamma, k) for m in m_values])
        orders = np.log(errors[:-1] / errors[1:]) / np.log(h_values[:-1] / h_values[1:])
        all_errors.append(errors)
        all_orders.append(orders)

    plot_convergence(k_values, h_values, all_errors, all_orders)
