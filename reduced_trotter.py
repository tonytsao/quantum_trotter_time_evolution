from qiskit import QuantumCircuit, Aer, execute
from qiskit.quantum_info.operators import Operator
import numpy as np
import os


from qmps.time_evolve_tools import put_env_on_right_site_D4, put_env_on_right_site
from qmps.ground_state import Hamiltonian
from scipy.linalg import expm
import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft


class F_cirq:
    def __init__(self, N, cnt):
        self.N = N
        self.cnt = cnt
        self.cb_idx = self.N

        def InitCirc():
            r_qu_n = self.N
            one_lay_gate_num = (2*self.N -1) // 2 ##one layer gate number
            #anci_qu_n = one_lay_gate_num * (cnt + 1)
            anci_qu_n = 0
            total_qubit = r_qu_n + anci_qu_n
            circ = QuantumCircuit(total_qubit)
            return circ

    
        self.__circ = InitCirc()

    def is_odd(self):
        return True if (self.N)%2 == 1 else False
    
    def Append_U(self, mat, idx):
        self.__circ.unitary(mat, idx)

    def Append_Z(self, idx):
        self.__circ.z(idx)

    def Append_U4(self, alpha, gamma, idx1, idx2, method = '2CNOT'):
        if method == '2CNOT':
            self.__circ.cx(idx1, idx2)
            if alpha != 0:
                self.__circ.rx(-2*alpha, idx1)
            self.__circ.rz(-2*gamma, idx2)
            self.__circ.cx(idx1, idx2)
        elif method == 'Rzx':
            self.__circ.ry(-(np.pi/2), idx2)
            self.__circ.rzx(2*gamma, idx1, idx2)
            self.__circ.ry(np.pi/2, idx2)


    def Append1Layer_U4(self, J, hx, hz, dt, mode, method):
        
        if mode == 1:
            for i in range(self.N):
                self.__circ.rx(-dt*hx, i)
                self.__circ.rz(-2*dt*hz, i)
        elif mode == 0:
            for i in range(1, 2):
                self.__circ.rz(-2*dt*hz, i)
                self.__circ.rx(-2*dt*hx, i)
                

        def Gate1List():
            gateN = self.N // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
    
        def Gate2List():
            gateN1 = self.N // 2
            gateN = (self.N - 1) // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i + 1
                #qa_idx = self.N + gateN1 + i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
        gate1list = Gate1List()
        gate2list = Gate2List()
        for gate_idxs in gate1list:
            #print('gate_idxs = ', gate_idxs)
            self.Append_U4(0, dt*J, gate_idxs[0], gate_idxs[1], method)

        for gate_idxs in gate2list:
            #print('gate_idxs = ', gate_idxs)
            self.Append_U4(0, dt*J, gate_idxs[0], gate_idxs[1], method)

        if mode == 1:
            for i in range(self.N):
                self.__circ.rx(-dt*hx, i)

    def Append1Layer_U4_dagger(self, J, hx, hz, dt, mode, method):
        
        if mode == 1:
            for i in range(self.N):
                self.__circ.rx(-dt*hx, i)
                #self.__circ.rz(-2*dt*hz, i)

        def Gate1List():
            gateN = self.N // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
    
        def Gate2List():
            gateN1 = self.N // 2
            gateN = (self.N - 1) // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i + 1
                #qa_idx = self.N + gateN1 + i
                lst.append([q1_idx, q1_idx + 1])
            #print(lst)
            return lst
        gate1list = Gate1List()
        gate2list = Gate2List()
        for gate_idxs in gate2list:
            #print('gate_idxs = ', gate_idxs)
            self.Append_U4(0, dt*J, gate_idxs[0], gate_idxs[1], method)

        for gate_idxs in gate1list:
            #print('gate_idxs = ', gate_idxs)
            self.Append_U4(0, dt*J, gate_idxs[0], gate_idxs[1], method)

        if mode == 1:
            for i in range(self.N):
                self.__circ.rz(-2*dt*hz, i)
                self.__circ.rx(-dt*hx, i)
        elif mode == 0:
            for i in range(1, 2):
                self.__circ.rx(-2*dt*hx, i)
                self.__circ.rz(-2*dt*hz, i)
                  
            


    def StateVectorRes(self, cnt):
        simulator = Aer.get_backend('statevector_simulator')
        #simulator = Aer.get_backend('qasm_simulator')
        result = execute(self.__circ, simulator).result()
        phi = result.get_statevector(self.__circ)
        #phi = phi/np.linalg.norm(phi)
        #print(phi)
        res = phi[0] #* pow(2, 6) / self.N
        return res

    def MeasRes2(self):

        #provider = IBMQ.get_provider(hub = 'ibm-q')
        #device = provider.get_backend('ibm_cairo')
        #simulator = Aer.get_backend('statevector_simulator')
        service = QiskitRuntimeService()
        backend = service.backend("ibm_osaka")
        options = Options()
        options.resilience_level = 1
        options.optimization_level = 3
        sampler = Sampler(backend=backend, options = options)
        shot_num = 4096
        job = sampler.run(self.__circ, shots = shot_num, initial_layout = [66, 67, 68, 69, 70])
        print(f">>> Job ID: {job.job_id()}")
        print(f">>> Job Status: {job.status()}")
        result = job.result()
        target = ''
        for i in range(self.cb_idx):
            target = target + '0'
        print(result)
        counts = result.quasi_dists[0][0]
        print(counts)
        #print(type(counts))
        #print(result.get_counts(self.__circ))
        #rate = np.sqrt(counts)
        return counts

    def PrintCirc(self):
        print(self.__circ)
    
    def PrintGates(self):
        print(dict(self.__circ.count_ops()))

def fqte_trot(N, J, hx, hz, cnt, cut, mode, method = '2CNOT'):
    dt = 0.1
    t = dt*cnt/cut
    #print(0.5*H)
    #tmp = put_env_on_left_site_D4(0.5*H)
    #tmp = tmp.reshape(2,2,2,2,2,2)
    #print(tmp[:,:,0,:,:,0].reshape(4,4)-H*0.5)
    #tmp = np.linalg.eig(H)
    #print(tmp[0])
    f_cirq = F_cirq(N, cnt)
    
    #print(H)
    #e_H = expm(-H*1j*dt)
    #print(e_H)
    #e_HU = Operator(e_H)#直接放
    #print("e_HU = ", e_HU)
    #e_H_bar = expm(H*1j*dt)
    #e_HU_bar = Operator(e_H_bar)
    #print("e_HU_bar = ", e_HU_bar)
    #tmp = put_env_on_right_site_D4(0.5*e_H)
    #tmp = tmp.reshape(2,2,2,2,2,2)
    #print(tmp)

    #inital state
    sigma_x = [[0, 1], [1, 0]]
    '''
    f_cirq.Append_U(sigma_x, 0)
    f_cirq.Append_U(sigma_x, 1)
    f_cirq.Append_U(sigma_x, 2)
    f_cirq.Append_U(sigma_x, 3)
    f_cirq.Append_U(sigma_x, 4)
    f_cirq.Append_U(sigma_x, 5)
    f_cirq.Append_U(sigma_x, 6)
    f_cirq.Append_U(sigma_x, 7)
    f_cirq.Append_U(sigma_x, 8)
    f_cirq.Append_U(sigma_x, 9)
    f_cirq.Append_U(sigma_x, 10)'''
    
    for i in range(0, cut):
        f_cirq.Append1Layer_U4(J, hx, hz, t, mode, method)

    sigma_z = [[1, 0], [0, -1]]

    mid_idx = N//2
    f_cirq.Append_Z(mid_idx)
    #f_cirq.Append_U(sigma_z, mid_idx)
    
    for i in range(0, cut):
        f_cirq.Append1Layer_U4_dagger(J, hx, hz, -t, mode, method)
    '''
    f_cirq.Append_U(sigma_x, 0)
    f_cirq.Append_U(sigma_x, 1)
    f_cirq.Append_U(sigma_x, 2)
    f_cirq.Append_U(sigma_x, 3)
    f_cirq.Append_U(sigma_x, 4)
    f_cirq.Append_U(sigma_x, 5)
    f_cirq.Append_U(sigma_x, 6)
    f_cirq.Append_U(sigma_x, 7)
    f_cirq.Append_U(sigma_x, 8)
    f_cirq.Append_U(sigma_x, 9)
    f_cirq.Append_U(sigma_x, 10)'''
    
    f_cirq.PrintCirc()
    f_cirq.PrintGates()
    
    val = f_cirq.StateVectorRes(cnt)
    #print(val)
    
    #f_cirq.PrintCirc()
    fac = pow(2, N -1)
    #e = (fac*val.real)
    #e = val.real
    er = val.real
    ei = val.imag
    e = (pow(er, 2)+pow(ei, 2))
    #print('magnetization', e)
    return e
    

## main function
#data_path = '/home/tonytsao/tfim/qte/'
def main(sampler=None):
    #g = 0.3
    #h =3
    #H0 = Hamiltonian({'ZZ':-1.0/2, 'X':g, 'Z':h})
    N = 3
    cnt = 50
    mode = 0
    method = '2CNOT'
    trotter_step = 7
    dt = 0.1
    es = []
    t = []
    J = -0.5
    hx = 0.3
    hz = 3
    f_e = open("t_simplify_7.txt", "w")
    for i in range(0, cnt):
        e = fqte_trot(N, J, hx, hz, i, trotter_step, mode, method)
        f_e.write(str(e) + '\n')
        f_e.flush()
        es.append(e)
        t.append(dt*i)
    #print(es)
    '''
    es1 = []
    for i in range(0, cnt):
        e = fqte_trot(N, J, hx, hz, i, trotter_step, mode, 'Rzx')
        #f_e.write(str(e) + '\n')
        #f_e.flush()
        es1.append(e)
        #t.append(dt*i)

    
    plt.plot(t, es,"-")
    plt.plot(t, es1,".")
    plt.title('Finite qte, initial state: DDDDDDD')
    plt.xlabel('t')
    plt.ylabel('sigma z')
    plt.savefig("2CNOT_rzx.jpg")
    plt.show()'''
    


if __name__ == "__main__":
    main()

