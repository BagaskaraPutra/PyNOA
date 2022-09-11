from sympy.algebras.quaternion import Quaternion

def quatsub_sym(q_sym_var, quat_sympy_num):
    # Function to convert sympy Quaternion object to sympy dictionary
    # Input args:
    # q_sym_var = symbolic quaternion in sympy vector
    # quat_sympy_num = numeric quaternion in sympy Quaternion object
    q_dict_num = {q_sym_var[0]: quat_sympy_num.a, 
                    q_sym_var[1]: quat_sympy_num.b, 
                    q_sym_var[2]: quat_sympy_num.c, 
                    q_sym_var[3]: quat_sympy_num.d}
    # print(q_dict_num)
    return q_dict_num