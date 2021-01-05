#Module Mapper for 2019 layout
def str_to_hw(string):
    string_to_hw = {}

    string_to_hw["L0T_axial"     ] = "F0H0"
    string_to_hw["L0T_stereo"    ] = "F0H1"
    string_to_hw["L1T_axial"     ] = "F0H2"
    string_to_hw["L1T_stereo"    ] = "F0H3"
    string_to_hw["L1B_axial"     ] = "F1H0"
    string_to_hw["L1B_stereo"    ] = "F1H1"
    string_to_hw["L0B_axial"     ] = "F1H2"
    string_to_hw["L0B_stereo"    ] = "F1H3"
    string_to_hw["L2T_axial"     ] = "F2H0"
    string_to_hw["L2T_stereo"    ] = "F2H1"
    string_to_hw["L3T_stereo"    ] = "F2H2"
    string_to_hw["L3T_axial"     ] = "F2H3"
    string_to_hw["L2B_stereo"    ] = "F3H2" #(real val F3H0)
    string_to_hw["L2B_axial"     ] = "F3H3" #swapped with 3B_axial (Real Val F3H1)
    string_to_hw["L3B_stereo"    ] = "F3H0" #(real val F3H0)
    string_to_hw["L3B_axial"     ] = "F3H1" #(real val F3H3)

   # string_to_hw["L2B_stereo"    ] = "F3H0" 
   # string_to_hw["L2B_axial"     ] = "F3H1" 
   # string_to_hw["L3B_stereo"    ] = "F3H2" 
   # string_to_hw["L3B_axial"     ] = "F3H3" 
    string_to_hw["L4T_axial_ele" ] = "F4H0"
    string_to_hw["L4T_axial_pos" ] = "F4H1"
    string_to_hw["L4T_stereo_ele"] = "F4H2"
    string_to_hw["L4T_stereo_pos"] = "F4H3"
    string_to_hw["L4B_stereo_ele"] = "F5H0"
    string_to_hw["L4B_stereo_pos"] = "F5H1"
    string_to_hw["L4B_axial_ele" ] = "F5H2"
    string_to_hw["L4B_axial_pos" ] = "F5H3"
    string_to_hw["L5T_axial_ele" ] = "F6H0"
    string_to_hw["L5T_axial_pos" ] = "F6H1"
    string_to_hw["L5T_stereo_ele"] = "F6H2"
    string_to_hw["L5T_stereo_pos"] = "F6H3"
    string_to_hw["L5B_stereo_ele"] = "F7H0"
    string_to_hw["L5B_stereo_pos"] = "F7H1"
    string_to_hw["L5B_axial_ele" ] = "F7H2"
    string_to_hw["L5B_axial_pos" ] = "F7H3"
    string_to_hw["L6T_axial_ele" ] = "F8H0"
    string_to_hw["L6T_axial_pos" ] = "F8H1"
    string_to_hw["L6T_stereo_ele"] = "F8H2"
    string_to_hw["L6T_stereo_pos"] = "F8H3"
    string_to_hw["L6B_stereo_ele"] = "F9H0"
    string_to_hw["L6B_stereo_pos"] = "F9H1"
    string_to_hw["L6B_axial_ele" ] = "F9H2"
    string_to_hw["L6B_axial_pos" ] = "F9H3"
    return string_to_hw[string]

def str_to_sw(string):
    string_to_sw = {}
    string_to_sw["L0T_axial"     ] = "ly1_m0" 
    string_to_sw["L0T_stereo"    ] = "ly2_m0" 
    string_to_sw["L1T_axial"     ] = "ly3_m0" 
    string_to_sw["L1T_stereo"    ] = "ly4_m0" 
    string_to_sw["L1B_axial"     ] = "ly4_m1" 
    string_to_sw["L1B_stereo"    ] = "ly3_m1" 
    string_to_sw["L0B_axial"     ] = "ly2_m1" 
    string_to_sw["L0B_stereo"    ] = "ly1_m1" 
    string_to_sw["L2T_axial"     ] = "ly5_m0" 
    string_to_sw["L2T_stereo"    ] = "ly6_m0" 
    string_to_sw["L3T_stereo"    ] = "ly8_m0" 
    string_to_sw["L3T_axial"     ] = "ly7_m0" 
    string_to_sw["L2B_stereo"    ] = "ly5_m1" 
    string_to_sw["L2B_axial"     ] = "ly6_m1" 
    string_to_sw["L3B_stereo"    ] = "ly7_m1" 
    string_to_sw["L3B_axial"     ] = "ly8_m1" 
    string_to_sw["L4T_axial_ele" ] = "ly9_m0" 
    string_to_sw["L4T_axial_pos" ] = "ly9_m2" 
    string_to_sw["L4T_stereo_ele"] = "ly10_m0"
    string_to_sw["L4T_stereo_pos"] = "ly10_m2"
    string_to_sw["L4B_stereo_ele"] = "ly9_m1" 
    string_to_sw["L4B_stereo_pos"] = "ly9_m3" 
    string_to_sw["L4B_axial_ele" ] = "ly10_m1"
    string_to_sw["L4B_axial_pos" ] = "ly10_m3"
    string_to_sw["L5T_axial_ele" ] = "ly11_m0"
    string_to_sw["L5T_axial_pos" ] = "ly11_m2"
    string_to_sw["L5T_stereo_ele"] = "ly12_m0"
    string_to_sw["L5T_stereo_pos"] = "ly12_m2"
    string_to_sw["L5B_stereo_ele"] = "ly11_m1"
    string_to_sw["L5B_stereo_pos"] = "ly11_m3"
    string_to_sw["L5B_axial_ele" ] = "ly12_m1"
    string_to_sw["L5B_axial_pos" ] = "ly12_m3"
    string_to_sw["L6T_axial_ele" ] = "ly13_m0"
    string_to_sw["L6T_axial_pos" ] = "ly13_m2"
    string_to_sw["L6T_stereo_ele"] = "ly14_m0"
    string_to_sw["L6T_stereo_pos"] = "ly14_m2"
    string_to_sw["L6B_stereo_ele"] = "ly13_m1"
    string_to_sw["L6B_stereo_pos"] = "ly13_m3"
    string_to_sw["L6B_axial_ele" ] = "ly14_m1"
    string_to_sw["L6B_axial_pos" ] = "ly14_m3"
    return string_to_sw[string]

def hw_to_str(string):
    hw_to_string = {}
    hw_to_string["F0H0"] = "L0T_axial"     
    hw_to_string["F0H1"] = "L0T_stereo"    
    hw_to_string["F0H2"] = "L1T_axial"     
    hw_to_string["F0H3"] = "L1T_stereo"    
    hw_to_string["F1H0"] = "L1B_axial"    
    hw_to_string["F1H1"] = "L1B_stereo"    
    hw_to_string["F1H2"] = "L0B_axial"     
    hw_to_string["F1H3"] = "L0B_stereo"    
    hw_to_string["F2H0"] = "L2T_axial"     
    hw_to_string["F2H1"] = "L2T_stereo"    
    hw_to_string["F2H2"] = "L3T_stereo"    
    hw_to_string["F2H3"] = "L3T_axial"     
    hw_to_string["F3H0"] = "L2B_stereo"    
    hw_to_string["F3H1"] = "L2B_axial"     
    hw_to_string["F3H2"] = "L3B_stereo"    
    hw_to_string["F3H3"] = "L3B_axial"     
    hw_to_string["F4H0"] = "L4T_axial_ele" 
    hw_to_string["F4H1"] = "L4T_axial_pos" 
    hw_to_string["F4H2"] = "L4T_stereo_ele"
    hw_to_string["F4H3"] = "L4T_stereo_pos"
    hw_to_string["F5H0"] = "L4B_stereo_ele"
    hw_to_string["F5H1"] = "L4B_stereo_pos"
    hw_to_string["F5H2"] = "L4B_axial_ele" 
    hw_to_string["F5H3"] = "L4B_axial_pos" 
    hw_to_string["F6H0"] = "L5T_axial_ele" 
    hw_to_string["F6H1"] = "L5T_axial_pos" 
    hw_to_string["F6H2"] = "L5T_stereo_ele"
    hw_to_string["F6H3"] = "L5T_stereo_pos"
    hw_to_string["F7H0"] = "L5B_stereo_ele"
    hw_to_string["F7H1"] = "L5B_stereo_pos"
    hw_to_string["F7H2"] = "L5B_axial_ele" 
    hw_to_string["F7H3"] = "L5B_axial_pos" 
    hw_to_string["F8H0"] = "L6T_axial_ele" 
    hw_to_string["F8H1"] = "L6T_axial_pos" 
    hw_to_string["F8H2"] = "L6T_stereo_ele"
    hw_to_string["F8H3"] = "L6T_stereo_pos"
    hw_to_string["F9H0"] = "L6B_stereo_ele"
    hw_to_string["F9H1"] = "L6B_stereo_pos"
    hw_to_string["F9H2"] = "L6B_axial_ele" 
    hw_to_string["F9H3"] = "L6B_axial_pos" 
    return hw_to_string[string]

