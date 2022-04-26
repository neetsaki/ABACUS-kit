import os
import numpy as np
def read_file(filepath):
    with open(filepath, 'r') as f:
        return f.readlines()
def KPT_gen(Gen_type,KPT_type,a,b,c,d_KPT,N=0,K_coord=[]):
    f_line1="K_POINTS\n"
    if Gen_type == 'auto':
        f_line2="0\n"
        a_bar=1/a
        b_bar=1/b
        c_bar=1/c
        N_a=int(max(a_bar/d_KPT,1))
        N_b=int(max(b_bar/d_KPT,1))
        N_c=int(max(c_bar/d_KPT,1))
        f_line4=str(N_a)+" "+str(N_b)+" "+str(N_c)+" 0 0 0" #total number of k-point, `0' means generate automatically
        if KPT_type == 'gamma':
            f_line3="Gamma\n"
        elif KPT_type == 'MP':
            f_line3="MP\n"
        else:
            print('KPT_type error')
    elif KPT_type == 'explicit':
        f_line2=str(N)+"\n"
        if KPT_type == 'Direct':
            f_line3='Direct\n'
        elif KPT_type == 'Cartesian':
            f_line3='Cartesian\n'
        for line_num in range(N):
            f_line4=str(K_coord[line_num][0])+" "+str(K_coord[line_num][1])+" "+str(K_coord[line_num][2])+"\n"
    else:
        print('KPT_type error')
    f=open("KPT","w")
    f.write(f_line1+f_line2+f_line3+f_line4)
    return True
#POSCAR可以借助pymatgen方便地建模，所以这里选用vasp的POSCAR格式做转化
def STRU_gen_from_POSCAR(POSCAR,PP_type="SG15_ONCV",method="pw",PP_path="/home/room/Software/QE-PP/SG15_ONCV_v1.0_upf",Orbital_path="/home/room/Software/QE-PP/Orb_SG15_DZP_E100_Standard_v1.0"):
    poscar_raw=read_file(POSCAR)
    atom_species_list=[]
    for atom in poscar_raw[5].strip ().split(" "):
        if atom == '':
            atom_species_list=atom_species_list
        else:
            atom_species_list.append(atom)
    atom_num_list=[]
    for atom_num in poscar_raw[6].strip ().split(" "):
        if atom_num == '':
            atom_num_list=atom_num_list
        else:
            atom_num_list.append(atom_num)
    s_line1="ATOMIC_SPECIES\n"
    #The mass of the elment is only used in molecular dynamics simulations.
    s_line1_1=""
    #The pseudopotential file is used to calculate the potential energy of the atoms.
    if method == "pw":
        for atom_species in atom_species_list:
            s_line1_1 += atom_species+"\t"+"1.0 "+atom_species+"_ONCV_PBE-1.0.upf"+"\n"
            PP_path="/home/room/Software/QE-PP/SG15_ONCV_v1.0_upf"
        s_line1_2="\n"
    elif method == "lcao":
#        Orbital_path="/home/room/Software/QE-PP/Orb_SG15_DZP_E100_Standard_v1.0"
        Orbital_list=os.listdir(Orbital_path)
        for atom_species in atom_species_list:
            s_line1_1 += atom_species+"\t"+"1.0 "+atom_species+"_ONCV_PBE-1.0.upf"+"\n"
        s_line1_2_1="NUMERICAL_ORBITAL"
        for atom_species in atom_species_list:
            for Orbital_file in Orbital_list:
                if atom_species == Orbital_file.split("_")[0]:
                    s_line1_2_1 += Orbital_file+"\n"
                else:
                    s_line1_2_1 = s_line1_2_1
        s_line1_2="\n"+s_line1_2_1+"\n"
    #this part ("lcao_reduced") is under developing. 
    elif method == "lcao_reduced":
        Orbital_path="/home/room/Software/QE-PP/Traditional_PP_Orb"
        Orbital_list=os.listdir(Orbital_path)
        for atom_species in atom_species_list:
            s_line1_1 += atom_species+"\t"+"1.0 "+atom_species+".PBE.UPF"+"\n"
        s_line1_2_1="NUMERICAL_ORBITAL"
        s_line1_2_2="\n"
    else:
        print("method error")
    s_line2="LATTICE_CONSTANT\n"
    s_line2_1=str(poscar_raw[1].strip ())+"\n"
    s_line3="LATTICE_VECTORS\n"
    s_line3_1=str(poscar_raw[2].strip ())+"\n"
    s_line3_2=str(poscar_raw[3].strip ())+"\n"
    s_line3_3=str(poscar_raw[4].strip ())+"\n\n"
    s_line4="ATOMIC_POSITIONS\n"
    if poscar_raw[7].strip()[0] == "C":
        s_line4_1="Cartesian\n"
        ini_line=8
    elif poscar_raw[7].strip()[0] == "D":
        s_line4_1="Direct\n"
        ini_line=8
    else:
        s_line4_1=poscar_raw[8].strip()
        ini_line=9 
        print("Check Fixed atoms !!")    #固定原子格式的转化尚未实现自动转化
    s_line4_2="\n"
    s_line5=""
    #coordinate start from 9th line

    for i in range(len(atom_num_list)):
        #如果需要改变原始磁矩，改变这里的0.0
        s_line5 += atom_species_list[i]+"\n"+"0.0\n"+str(atom_num_list[i])+"\n"
        #写入坐标，默认放开所有原子
        for j in range(int(atom_num_list[i])):
            s_line5 += str(poscar_raw[8+j].strip ())+" 1 1 1 \n"
        ini_line += int(atom_num_list[i])
        s_line5 += "\n"
    s_group_1=s_line1+s_line1_1+s_line1_2
    s_group_2=s_line2+s_line2_1
    s_group_3=s_line3+s_line3_1+s_line3_2+s_line3_3
    s_group_4=s_line4+s_line4_1+s_line4_2
    s_group_5=s_line5
    s=open("STRU","w")
    s.write(s_group_1+s_group_2+s_group_3+s_group_4+s_group_5)
    s.close()
    return atom_species_list,atom_num_list,s_line3_1,s_line3_2,s_line3_3

def INPUT_gen(model,PP_PATH,atom_type_num,N_band,smearing=True):
    i_line0="INPUT_PARAMETERS\n"
    if model == "scf_pw":
    #General parameters
        i_line1_0="calculation\t\t\t"+"scf"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"pw"+"\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line1_2+i_line1_3+i_line1_4+i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6
        input_line=i_group1+i_group2
        i=open("INPUT","w")
        i.write(input_line)
        i.close()
    elif model == "scf_lcao":
    #General parameters
        i_line1_0="calculation\t\t\t"+"scf"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"lcao"+"\n"+"gamma_only\t\t\t"+"1\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line1_2+i_line1_3+i_line1_4+i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6
        input_line=i_group1+i_group2
        i=open("INPUT","w")
        i.write(input_line)
        i.close()
    elif model == "relax_pw":
    #General parameters
        i_line1_0="calculation\t\t\t"+"relax"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"pw"+"\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6+"\n"
    #Geometry relax
        i_line3_1="relax_nmax\t\t\t"+"100\n"
        i_line3_2="force_thr_ev\t\t\t"+"1.0e-3\n"
        i_line3_3="relax_method\t\t\t"+"bfgs\n"
        i_line3_4="fixed_axes\t\t\t"+"None\n"
        i_line3_5="cal_stress\t\t\t"+"0\n"
        i_group3=i_line3_1+i_line3_2+i_line3_3+i_line3_4+i_line3_5
        input_line=i_group1+i_group2+i_group3
        i=open("INPUT","w")
        i.write(input_line)
        i.close()
    elif method == "relax_lcao":
    #General parameters
        i_line1_0="calculation\t\t\t"+"relax"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"lcao"+"\n"+"gamma_only\t\t\t"+"1\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6+"\n"
    #Geometry relax
        i_line3_1="relax_nmax\t\t\t"+"100\n"
        i_line3_2="force_thr_ev\t\t\t"+"1.0e-3\n"
        i_line3_3="relax_method\t\t\t"+"bfgs\n"
        i_line3_4="fixed_axes\t\t\t"+"None\n"
        i_line3_5="cal_stress\t\t\t"+"0\n"
        i_group3=i_line3_1+i_line3_2+i_line3_3+i_line3_4+i_line3_5
        input_line=i_group1+i_group2+i_group3
        i=open("INPUT","w")
        i.write(input_line)
        i.close()
    elif method == "dos_lcao":
    #General parameters
        i_line1_0="calculation\t\t\t"+"relax"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"lcao"+"\n"+"gamma_only\t\t\t"+"1\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line1_2+i_line1_3+i_line1_4+i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6+"\n"
    #Geometry relax
        i_line3_1="relax_nmax\t\t\t"+"100\n"
        i_line3_2="force_thr_ev\t\t\t"+"1.0e-3\n"
        i_line3_3="relax_method\t\t\t"+"bfgs\n"
        i_line3_4="fixed_axes\t\t\t"+"None\n"
        i_line3_5="cal_stress\t\t\t"+"0\n"
        i_group3=i_line3_1+i_line3_2+i_line3_3+i_line3_4+i_line3_5
        i_line4="out_chg\t\t\t"+"1\n"
        input_line=i_group1+i_group2+i_group3+i_line4
        #输出第一步结构优化的INPUT文件
        i=open("INPUT1-relax","w")
        i.write(input_line)
        i.close()
        i_line1_0="calculation\t\t\t"+"nscf"+"\n"
        i_line4_1="init_chg\t\t\t"+"file\n"
        i_line4_2="out_dos\t\t\t"+"1\n"
        i_group4=i_line4_1+i_line4_2
        input_line=i_group1+i_group2+i_group3+i_group4
        #输出第二步DOS的INPUT文件
        i=open("INPUT2-DOS","w")
        i.write(input_line)
        i.close()
    elif method == "dos_pw":
    #General parameters
        i_line1_0="calculation\t\t\t"+"relax"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"pw"+"\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line1_2+i_line1_3+i_line1_4+i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6+"\n"
    #Geometry relax
        i_line3_1="relax_nmax\t\t\t"+"100\n"
        i_line3_2="force_thr_ev\t\t\t"+"1.0e-3\n"
        i_line3_3="relax_method\t\t\t"+"bfgs\n"
        i_line3_4="fixed_axes\t\t\t"+"None\n"
        i_line3_5="cal_stress\t\t\t"+"0\n"
        i_group3=i_line3_1+i_line3_2+i_line3_3+i_line3_4+i_line3_5
        i_line4="out_chg\t\t\t"+"1\n"
        input_line=i_group1+i_group2+i_group3+i_line4
        #输出第一步结构优化的INPUT文件
        i=open("INPUT1-relax","w")
        i.write(input_line)
        i.close()
        i_line1_0="calculation\t\t\t"+"nscf"+"\n"
        i_line4_1="init_chg\t\t\t"+"file\n"
        i_line4_2="out_dos\t\t\t"+"1\n"
        i_group4=i_line4_1+i_line4_2
        input_line=i_group1+i_group2+i_group3+i_group4
        #输出第二步DOS的INPUT文件
        i=open("INPUT2-DOS","w")
        i.write(input_line)
        i.close()
    elif method == "band_pw":
    #General parameters
        i_line1_0="calculation\t\t\t"+"relax"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"pw"+"\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line1_2+i_line1_3+i_line1_4+i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6
    #Geometry relax
        i_line3_1="relax_nmax\t\t\t"+"100\n"
        i_line3_2="force_thr_ev\t\t\t"+"1.0e-3\n"
        i_line3_3="relax_method\t\t\t"+"bfgs\n"
        i_line3_4="fixed_axes\t\t\t"+"None\n"
        i_line3_5="cal_stress\t\t\t"+"0\n"
        i_group3=i_line3_1+i_line3_2+i_line3_3+i_line3_4+i_line3_5
        i_line4="out_chg\t\t\t"+"1\n"
        input_line=i_group1+i_group2+i_group3+i_line4
        #输出第一步结构优化的INPUT文件
        i=open("INPUT1-relax","w")
        i.write(input_line)
        i.close()
        i_line1_0="calculation\t\t\t"+"nscf"+"\n"
        i_line4_1="init_chg\t\t\t"+"file\n"
        i_line4_2="out_band\t\t\t"+"1\n"
        i_group4=i_line4_1+i_line4_2
        input_line=i_group1+i_group2+i_group3+i_group4
        #输出第二步BAND的INPUT文件
        i=open("INPUT2-BAND","w")
        i.write(input_line)
        i.close()
    elif method == "band_lcao":
    #General parameters
        i_line1_0="calculation\t\t\t"+"relax"+"\n"
        i_line1_1="pseudo_dir\t\t\t"+PP_PATH+"\n"
        i_line1_2="ntype\t\t\t"+str(atom_type_num)+"\n"
        i_line1_3="nbands\t\t\t"+str(N_band)+"\n"
        i_line1_4="basis_type\t\t\t"+"lcao"+"\n"+"gamma_only\t\t\t"+"1\n"
        i_group1=i_line0+i_line1_0+i_line1_1+i_line1_2+i_line1_3+i_line1_4+"\n"
    #Self-consistent calculation
        i_line2_1="ecutwfc\t\t\t"+"50.0"+"\n"
        i_line2_2="scf_nmax\t\t\t"+"60"+"\n"
        i_line2_3="symmetry\t\t\t"+"0"+"\n"
        i_line2_4="scf_thr\t\t\t"+"1.0E-8"+"\n"
        if smearing==True:
            i_line2_5="smearing_method\t\t"+"Gaussian"+"\n"
            i_line2_6="smearing_sigma\t\t"+"0.021"+"\n"
        else:
            i_line2_5=""
            i_line2_6=""
        i_group2=i_line1_2+i_line1_3+i_line1_4+i_line2_1+i_line2_2+i_line2_3+i_line2_4+i_line2_5+i_line2_6+"\n"
    #Geometry relax
        i_line3_1="relax_nmax\t\t\t"+"100\n"
        i_line3_2="force_thr_ev\t\t\t"+"1.0e-3\n"
        i_line3_3="relax_method\t\t\t"+"bfgs\n"
        i_line3_4="fixed_axes\t\t\t"+"None\n"
        i_line3_5="cal_stress\t\t\t"+"0\n"
        i_group3=i_line3_1+i_line3_2+i_line3_3+i_line3_4+i_line3_5
        i_line4="out_chg\t\t\t"+"1\n"
        input_line=i_group1+i_group2+i_group3+i_line4
        #输出第一步结构优化的INPUT文件
        i=open("INPUT1-relax","w")
        i.write(input_line)
        i.close()
        i_line1_0="calculation\t\t\t"+"nscf"+"\n"
        i_line4_1="init_chg\t\t\t"+"file\n"
        i_line4_2="out_dos\t\t\t"+"1\n"
        i_group4=i_line4_1+i_line4_2
        input_line=i_group1+i_group2+i_group3+i_group4
        #输出第二步BAND的INPUT文件
        i=open("INPUT2-BAND","w")
        i.write(input_line)
        i.close()
    print("input file generated!")
    return True

PP_path_custom="./SG15_ONCV_v1.0_upf"
BS_path_custom="./Orb_SG15_DZP_E60_Full_v1.0/Orb_SG15_DZP_E60_Full_v1.0"
INPUT_file_Template="./INPUT_Template"

[atom_spec,atom_nu,a,b,c]=STRU_gen_from_POSCAR("./biphen-P.vasp","SG15_ONCV","pw",PP_path_custom,BS_path_custom)
a=[float(a_i) for a_i in a.strip().split("         ")]
Lattice_a=np.sqrt(a[0]**2+a[1]**2+a[2]**2)
b=[float(b_i) for b_i in b.strip().split("        ")]
Lattice_b=np.sqrt(b[0]**2+b[1]**2+b[2]**2)
c=[float(c_i) for c_i in c.strip().split("        ")]
Lattice_c=np.sqrt(c[0]**2+c[1]**2+c[2]**2)
KPT_gen('auto','gamma',Lattice_a,Lattice_b,Lattice_c,0.04)
total_atom_num=sum([float(nu) for nu in atom_nu])
INPUT_gen("relax_pw",PP_path_custom,len(atom_spec),int(8*total_atom_num))
