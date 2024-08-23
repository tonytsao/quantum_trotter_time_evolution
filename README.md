Quantum Trotter Time Evolution
===
In the repository, I store quantum trotter time evolution code. The physics model we simulated is Transverse Field Ising Model.  
**trotter.py**: Quantum simulation using Trotter-Suzuki decomposition with fixed time step.
**reduced_trotter.py**: In order to decrease the needed two-qubits gate, we use reduced Trotter block instead.
**real_device_reduced_trotter.py**: Put reduced_trotter.py on IBMQ device.
