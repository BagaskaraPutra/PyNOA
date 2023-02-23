import sympy as sp
from sympy.solvers import pdsolve, checkpdesol
import random
import json
from itertools import permutations, combinations, product
from IPython.display import display

from datetime import datetime
import os
import pickle, dill

class NOA():
    def __init__(self, name):
        self.name = name                        # object name
        self.x = sp.Matrix()                    # state vector
        self.f = []                             # vector fields
        self.h = sp.Matrix()                    # measurement model
        self.sys_order = int(1)                 # system order / number of states
        self.num_inputs = int(1)                # number of control inputs
        
        self.json_config_name = ""              # json parameters configuration file name
        self.params_config_subs = sp.Matrix()   # parameters that will be substituted with json numerical values
        self.new_params_dict = {}

        self.combn_permn_opt = "permutation"    # combination or permutation of vector fields Lie derivatives
        self.LD_order = int(0)                  # maximum Lie derivative order
        self.Lfh = {}                           # Lie derivatives
        self.dLfh_dx = {}                       # gradient of Lie derivatives wrt x
        
        self.rank_calc_opt = "symbolic"         # symbolic or numeric rank calculation of observability matrix
        self.numeric_params_dict = {}           # numeric parameters dictionary
        self.obsv_mat = sp.Matrix()             # observability matrix
        self.obsv_mat_num = sp.Matrix()         # numeric observability matrix substituted with json numerical values
        self.rank_obsv_mat = int(0)             # rank of observability matrix

        self.null_calc_opt = "symbolic"         # symbolic or numeric nullspace calculation of observability matrix
        self.cont_symm = []                     # continuous symmetries (nullspace of observability matrix)
        self.cont_symm_mat = sp.Matrix()        # matrix comprised of continuous symmetries
        self.nonobsv_subspace = sp.Matrix()     # joint observable & unobservable subspace
        self.obsv_subspace = sp.Matrix()        # observable subspace

        self.backup_name = ""                   # backup pickle & dill file name
        print(name, "NOA object initialized")

    def permn_rep(self, input_list, r_length_permn):
        # Permutation with repetition method
        # input_list     : list of elements to be permutated
        # r_length_permn : r length permutations of elements
        results = []
        for prod in product(input_list, repeat = r_length_permn):
            results.append(prod)
        return results
    
    def update_params_dict(self, new_params_dict):
        self.new_params_dict = new_params_dict

    def obsv_mat_construct(self, idx_all_perm, k):
        for idx_perm_k_iter in idx_all_perm:
            idx_perm_k = [*idx_perm_k_iter]
            # print(idx_perm_k)
            previous_order_vector_field_str = "k" + str(k-1)
            current_order_vector_field_str  = "k" + str(k)
            
            for idx_vector_field in idx_perm_k[:-1]:
                if len(idx_perm_k) > 1:
                    previous_order_vector_field_str = previous_order_vector_field_str + "f" + str(idx_vector_field)
            
            for idx_vector_field in idx_perm_k: 
                current_order_vector_field_str  = current_order_vector_field_str + "f" + str(idx_vector_field)
            
            self.Lfh[current_order_vector_field_str]     = self.dLfh_dx[previous_order_vector_field_str] * self.f[idx_vector_field]
            self.dLfh_dx[current_order_vector_field_str] = self.Lfh[current_order_vector_field_str].jacobian(self.x)
            gradient_rows = self.dLfh_dx[current_order_vector_field_str].rows
            gradient_cols = self.dLfh_dx[current_order_vector_field_str].cols
            if(self.dLfh_dx[current_order_vector_field_str] == sp.zeros(gradient_rows,gradient_cols)):
                print("current Lie derivative: ", current_order_vector_field_str, " not appended to observability matrix due to null vector")
            else:
                self.obsv_mat = self.obsv_mat.col_join(self.dLfh_dx[current_order_vector_field_str])
                print("current Lie derivative: ", current_order_vector_field_str, "appended to observability matrix")

    def ORC(self):
        # Observability Rank Criterion (ORC) for Nonlinear Observability Analysis (NOA)
        #
        # Input args:
        # x = state variables [x1; x2; ... xN]
        # f = control-affine process model [f1; f2; ... fM]
        # h = measurement model [h1; h2; ... hP];
        #
        # Output args:
        # Lfh      = Lie derivatives
        # dLfh_dx  = gradient of Lie derivatives wrt x
        # obsv_mat = observability matrix
        #
        # Reference: Continuous Symmetries and Observability Properties in Autonomous Navigation (Martinelli, 2010)

        # Observability rank criterion
        self.sys_order = self.x.rows      # order of system / number of states
        self.num_inputs = len(self.f)-1   # number of inputs
        if (self.LD_order == 0):
            self.LD_order = self.num_inputs

        # zero-th order Lie derivative & its gradient wrt x
        self.Lfh = {"k0": self.h}
        self.dLfh_dx = {"k0": self.Lfh["k0"].jacobian(self.x)}
        self.obsv_mat = self.dLfh_dx["k0"]                # initialize observability matrix

        # k-th order Lie derivative & its gradient wrt x
        if self.combn_permn_opt == "combination":
            print("Combination of Vector Fields")
            for k in range(1,self.LD_order+1):
                idx_all_perm = [*combinations(list(range(self.num_inputs+1)), k)]
                print(idx_all_perm)
                self.obsv_mat_construct(idx_all_perm, k)
            while (self.obsv_mat.rows < self.obsv_mat.cols):
                self.LD_order += 1
                if (self.LD_order < self.num_inputs):
                    print("Insufficient obv_mat rows. Appending Lie derivative order ", self.LD_order, " to obsv_mat ...")
                    idx_all_perm = [*combinations(list(range(self.num_inputs+1)), self.LD_order)]
                    print(idx_all_perm)
                    self.obsv_mat_construct(idx_all_perm, self.LD_order)
                else:
                    break
            print("Insufficient obsv_mat rows. Auto-construct obsv_mat using Permutation of vector fields ...")
            self.combn_permn_opt = "permutation"
            self.ORC()  
                
        
        elif self.combn_permn_opt == "drift2ndOrder":
            print("Drift 2nd Order Vector Fields")
            for k in range(1, 2+1):
                if k==1:
                    idx_all_perm = [*combinations(list(range(self.num_inputs+1)), k)]
                else:
                    idx_all_perm = []
                    for i2 in range(1, self.num_inputs+1):
                        idx_all_perm.append((0, i2))
                print(idx_all_perm)
                self.obsv_mat_construct(idx_all_perm, k)
            if (self.obsv_mat.rows < self.obsv_mat.cols):  
                print("Insufficient obsv_mat rows. Auto-construct obsv_mat using Permutation of vector fields ...")
                self.combn_permn_opt = "permutation"
                self.ORC()          
        
        elif self.combn_permn_opt == "drift2ndOrderWuest":
            print("Drift 2nd Order Vector Fields Wuest 2019")
            for k in range(1, 2+1):
                if k==1:
                    idx_all_perm = [*combinations(list(range(self.num_inputs+1)), k)]
                else:
                    idx_all_perm = []
                    for i2 in range(1, self.num_inputs+1):
                        idx_all_perm.append((i2, 0))
                print(idx_all_perm)
                self.obsv_mat_construct(idx_all_perm, k)
            if (self.obsv_mat.rows < self.obsv_mat.cols):  
                print("Insufficient obsv_mat rows. Auto-construct obsv_mat using Permutation of vector fields ...")
                self.combn_permn_opt = "permutation"
                self.ORC()   
        
        else:
            print("Permutation of Vector Fields")
            for k in range(1,self.LD_order+1):
                idx_all_perm = [*self.permn_rep(list(range(self.num_inputs+1)), k)]
                print(idx_all_perm)
                self.obsv_mat_construct(idx_all_perm, k)
            while (self.obsv_mat.rows < self.obsv_mat.cols):
                self.LD_order += 1
                print("Insufficient obv_mat rows. Appending Lie derivative order ",self.LD_order, " to obsv_mat ...")
                idx_all_perm = [*self.permn_rep(list(range(self.num_inputs+1)), self.LD_order)]
                print(idx_all_perm)
                self.obsv_mat_construct(idx_all_perm, self.LD_order)
                
        if self.rank_calc_opt == "numeric" or self.rank_calc_opt == "numerical":
            if self.json_config_name != "":
                with open(self.json_config_name) as json_file:
                    data = json.load(json_file)
                # numeric_params_list = []
                
                # Get numeric_params from json config file
                for keywords_params in self.params_config_subs:
                    numeric_params_get = data.get(str(keywords_params))
                    if numeric_params_get != None:
                        numeric_params_tuple = (str(keywords_params), data[str(keywords_params)])
                        self.numeric_params_dict[keywords_params] = data[str(keywords_params)]
                    else:
                        self.numeric_params_dict[keywords_params] = 0.0
                        print(f"{str(keywords_params)} is not found in {self.json_config_name}")
                    # numeric_params_list.append(numeric_params_tuple)
                print("Substituted parameters with values from ", self.json_config_name)
            else:
                first_primes = []
                number_of_first_primes = self.params_config_subs.rows
                num = 1
                while(len(first_primes) < number_of_first_primes):
                    num += 1
                    if num > 1:
                        for i in range(2, num):
                            if (num % i) == 0:
                                break
                        else:
                            first_primes.append(num)
                random.shuffle(first_primes)
                
                num = 0
                for keywords_params in self.params_config_subs:
                    self.numeric_params_dict[keywords_params] = first_primes[num]
                    num += 1   
                print("Substituted parameters with values from same-order prime numbers") 

            # Update numeric_params_dict if there are new parameters      
            if (self.new_params_dict != {}):
                if self.numeric_params_dict != {}:
                    self.numeric_params_dict = {key: self.new_params_dict.get(key, self.numeric_params_dict[key]) for key in self.numeric_params_dict}
                else:
                    self.numeric_params_dict = {key: self.new_params_dict.get(key, self.new_params_dict[key]) for key in self.new_params_dict}
                print("Updated numeric parameters dictionary!")
            print("printing self.numeric_params_dict:")
            print(self.numeric_params_dict)
            # self.obsv_mat_num = self.obsv_mat.subs(numeric_params_list)
            self.obsv_mat_num = self.obsv_mat.xreplace(self.numeric_params_dict)
            print("\nCalculating rank of numerical observability matrix ...")
            self.rank_obsv_mat = self.obsv_mat_num.rank()    
        else:
            print("\nCalculating rank of symbolic observability matrix ...")
            self.rank_obsv_mat = self.obsv_mat.rank()
        
        if(self.rank_obsv_mat == self.sys_order):
            print("The system is weakly locally observable (WLO)")
        else:
            print("The system is NOT weakly locally observable (WLO)!"); 
            print("Dimension of largest WLO subsystem (Observable & Jointly Observable): ", self.rank_obsv_mat)
            print("Dimension of Undistinguishable Region (Unobservable): ", self.sys_order-self.rank_obsv_mat)            

    # Observable, joint observable, and unobservable states
    def observable_mode(self):
        if self.null_calc_opt == "numeric" or self.null_calc_opt == "numerical":
            self.cont_symm = self.obsv_mat_num.nullspace()  # numerical continous symmetries
        else:
            self.cont_symm = self.obsv_mat.nullspace()      # symbolic continous symmetries
        
        if len(self.cont_symm) == 0:
            print("All states are observable")
            display(self.x.T)
        else:
            print("The observable modes are g(x) which satisfiy these partial differential equations (PDE): ")
            dg_latex = sp.Symbol('{\partial g(\mathbf{x})}')
            par_diff_latex = sp.symbols('\partial')
            for ws in self.cont_symm:
                pde_latex = 0
                for idx_x in range(self.sys_order):
                    pde_latex += ws[idx_x]*dg_latex/(par_diff_latex*self.x[idx_x])
                display(sp.Eq(pde_latex, 0))

            all_standard_basis      = True
            list_of_list_zero_idx   = [[]]*len(self.cont_symm)
            list_of_set_zero_idx    = [[]]*len(self.cont_symm)
            set_obsv_idx            = []

            # Check if the continuous symmetries are all standard basis to reduce computation
            for idx_ws, ws in enumerate(self.cont_symm):
                zero_counter = int(0)
                one_counter = int(0)
                list_of_list_zero_idx[idx_ws] = []
                
                for i, ws_i in enumerate(ws):
                    if ws_i == 0:
                        zero_counter += 1
                        # Mark which elements are 0 and compare them to other cont_symms to get observable modes    
                        list_of_list_zero_idx[idx_ws].append(i)
                    else:
                        one_counter += 1
                
                list_of_set_zero_idx[idx_ws] = set(list_of_list_zero_idx[idx_ws])
                
                if (one_counter == 1 and zero_counter == self.sys_order-1):
                    all_standard_basis = all_standard_basis and True
                else:
                    all_standard_basis = all_standard_basis and False

            for idx_ws, _ in enumerate(list_of_set_zero_idx):
                if idx_ws == 0:
                    set_obsv_idx = list_of_set_zero_idx[idx_ws]
                    if len(list_of_set_zero_idx) == 1:
                        break
                if idx_ws < len(list_of_set_zero_idx)-1:
                    set_obsv_idx = set_obsv_idx.intersection(list_of_set_zero_idx[idx_ws+1])
            
            # print("list_of_set_zero_idx: ", list_of_set_zero_idx)
            # print("set_obsv_idx: ", set_obsv_idx)

            if(all_standard_basis):
                self.cont_symm_mat      = sp.Matrix([self.cont_symm])
                self.nonobsv_subspace   = self.cont_symm_mat.T*self.x  # unobservable states
                set_nonobsv_subspace    = set(self.nonobsv_subspace)
                self.obsv_subspace      = set(self.x) - set_nonobsv_subspace
                print("Observable states: ")
                display(self.obsv_subspace)
                print("Unobservable states: ")
                display(set_nonobsv_subspace)   
            else:
                for ws in self.cont_symm:
                    g = sp.Function('g')
                    eval_str = 'g('
                    for idx_x in range(self.sys_order):
                        eval_str += 'self.x[' + str(idx_x) + ']'
                        if idx_x == self.sys_order-1:
                            eval_str += ')'
                            break
                        else:
                            eval_str += ','
                    u = eval(eval_str)
                    # ux = u.diff(x)
                    # uy = u.diff(y)
                    # genform = a*ux + b*uy + c*u
                    genform = 0
                    for idx_x in range(self.sys_order):
                        genform += ws[idx_x]*u.diff(self.x[idx_x])
                    # TODO: Solve partial differential equation for more than 2 state variables in sympy
                    # sol = pdsolve(genform)
                    print("Observable modes: ")
                    for obsv_idx in set_obsv_idx:
                        sol = self.x[obsv_idx]
                        check_sol = (checkpdesol(genform, sol))
                        if check_sol[0] == True:
                            display(sol)
    
    def save(self, dir):
        now = datetime.now()
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.backup_name = dir + "/" + self.name +"_"+ str(self.LD_order) + "OrderLD_" + self.combn_permn_opt + "_" + str(now.strftime("%Y-%m-%d_%H-%M-%S"))
        with open(self.backup_name+".pkl", 'wb') as file:
            pickle.dump(self, file)
        print("Pickle file saved to " + self.backup_name + ".pkl")
        dill_dumps = dill.dumps(self)
        dill.dump_session(self.backup_name + ".db")
        print("Dill file saved to " + self.backup_name + ".db")