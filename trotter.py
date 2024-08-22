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
    def __init__(self, N, H, cnt):
        self.N = N
        self.H = H
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

    def Append1Layer_U(self, mat):
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
            self.__circ.unitary(mat, gate_idxs)

        for gate_idxs in gate2list:
            #print('gate_idxs = ', gate_idxs)
            self.__circ.unitary(mat, gate_idxs)

    def Append1Layer_U_dagger(self, mat):
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
            self.__circ.unitary(mat, gate_idxs)

        for gate_idxs in gate1list:
            #print('gate_idxs = ', gate_idxs)
            self.__circ.unitary(mat, gate_idxs)


    def Append1Layer(self, mat):
        def Gate1List():
            gateN = self.N // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i
                lst.append([q1_idx, q1_idx + 1, self.cb_idx])
                self.cb_idx = self.cb_idx + 1
            #print(lst)
            return lst
    
        def Gate2List():
            gateN1 = self.N // 2
            gateN = (self.N - 1) // 2
            lst = []
            for i in range(0, gateN):
                q1_idx = 2 * i + 1
                qa_idx = self.N + gateN1 + i
                lst.append([q1_idx, q1_idx + 1, self.cb_idx])
                self.cb_idx = self.cb_idx + 1
            #print(lst)
            return lst
        gate1list = Gate1List()
        gate2list = Gate2List()
        for gate_idxs in gate1list:
            #print('gate_idxs = ', gate_idxs)
            self.__circ.unitary(mat, gate_idxs)

        for gate_idxs in gate2list:
            self.__circ.unitary(mat, gate_idxs)

        #self.__circ.reset(2)


    def StateVectorRes(self, cnt):
        simulator = Aer.get_backend('statevector_simulator')
        #simulator = Aer.get_backend('qasm_simulator')
        result = execute(self.__circ, simulator).result()
        phi = result.get_statevector(self.__circ)
        #phi = phi/np.linalg.norm(phi)
        #print(phi)
        res = phi[0] #* pow(2, 6) / self.N
        return res

    def PrintCirc(self):
        print(self.__circ)

def fqte(N, H, cnt, dt = 0.1):
    dt = dt
    #print(0.5*H)
    #tmp = put_env_on_left_site_D4(0.5*H)
    #tmp = tmp.reshape(2,2,2,2,2,2)
    #print(tmp[:,:,0,:,:,0].reshape(4,4)-H*0.5)
    #tmp = np.linalg.eig(H)
    #print(tmp[0])
    f_cirq = F_cirq(N, H, cnt)
    
    #print(H)
    e_H = expm(-H*1j*dt)
    #print(e_H)
    e_HU = Operator(e_H)#直接放
    #print("e_HU = ", e_HU)
    e_H_bar = expm(H*1j*dt)
    e_HU_bar = Operator(e_H_bar)
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
    
    for i in range(0, cnt):
        f_cirq.Append1Layer_U(e_HU)
    #f_cirq.Append1Layer(HU)

    sigma_z = [[1, 0], [0, -1]]

    mid_idx = N//2
    f_cirq.Append_U(sigma_z, mid_idx)
    #f_cirq.Append_U(sigma_x, mid_idx)
    
    for i in range(0, cnt):
        f_cirq.Append1Layer_U_dagger(e_HU_bar)
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
    
    #f_cirq.PrintCirc()
    
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

def trotter(N, J, hx, hz, cnt, dt = 0.1):
    g = hx
    h = hz
    H0 = Hamiltonian({'ZZ':J, 'X':g, 'Z':h})
    #print(H0.to_matrix())
    #dt = 0.1
    es = []
    t = []
    #f_e = open(data_path + "trotter.txt", "w")
    for i in range(0, cnt):
        e = fqte(N, H0.to_matrix(), i, dt)
        #f_e.write(str(e) + '\n')
        #f_e.flush()
        es.append(e)
        t.append(dt*i)
    #print(es)
    '''
    plt.plot(t, es,"-")
    plt.title('Finite qte, initial state: DDD')
    plt.xlabel('t')
    plt.ylabel('sigma z')
    plt.savefig("old_0.75_0.75.jpg")
    plt.show()'''
    return es

e = trotter(7, -0.5, 0.3, 3, 50, 0.1)
print(e)
'''
e1 = trotter(7, -0.25, 0.25, 3, 200)
print(e1)
t = []
t1 = []
for i in range(200):
    t.append(i*0.025)
    t1.append(i*0.025)


plt.plot(t, e, "-", label = '-1, 1, 12')
plt.plot(t1, e1, ".", label = '-0.25, 0.25, 3')
plt.title('Finite qte, initial state: DDDDDDD')
plt.xlabel('t')
plt.ylabel('$|\sigma_{z}|^{2}$')
plt.legend(loc = 'best')
plt.savefig("a.jpg")
plt.show()'''

#data_path = '/home/tonytsao/tfim/qte/data/'
f_e = open("trotter0.5.txt", "w")
for i in range(0, len(e)):
    f_e.write(str(e[i]) + '\n')
    f_e.flush()



