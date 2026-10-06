import plot_v2_node_decomposition
import plot_v2_scaling_by_observables
import plot_v2_speedup_decomposition
import plot_v2_strong_scaling
import plot_v2_threads
import plot_v2_weak_scaling

if __name__ == '__main__':
    for m in (plot_v2_threads, plot_v2_strong_scaling, plot_v2_scaling_by_observables,
              plot_v2_speedup_decomposition, plot_v2_node_decomposition, plot_v2_weak_scaling):
        m.main()
