"""
Strong Branching Criterion Simulation Script

This script simulates the strong branching process within a branch-and-bound tree to determine 
the optimal point for stopping strong branching, with the goal of minimizing the tree size.

Usage:
    - Configure the 'file_path' variable to specify the path to the input text file containing dual gain data.
    - Customize simulation parameters such as 'gap', 'num_repeats', and others as needed.

Description:
    The script performs the following steps:
    1. Reads variables and their dual gains from a text file.
    2. Assumes that dual gains follow an exponential distribution.
    3. Acquires an initial sample to estimate the parameters of this distribution.
    4. Computes the current tree size based on the variables selected thus far.
    5. Estimates the expected tree size if one additional iteration of strong branching is conducted.
    6. Determines the optimal stopping point for strong branching based on a predefined criterion.
"""

import math
import re 
import random
import numpy as np
from typing import List, Dict, Tuple

class Node:  
    def __init__(self, gap_closed):
        
        """
        Initialize a Node in a tree.

        Args:
            gap_closed (float): The gap value closed by this Node.
        """

        self.left = None
        self.right = None
        self.gap_closed = gap_closed
        
class Variable:

    """
    This class defines a basic representation for variables involved in strong branching.
    Variables can have left and right gains, which are used in branching decisions.

    Attributes:
        left_gain (float): The gain associated with branching left on this variable.
        right_gain (float): The gain associated with branching right on this variable.
    """

    def __init__(self):
        pass

    def get_gains(self) -> Tuple[float, float]:
        return (self.left_gain, self.right_gain)
    
    def __str__(self) -> str:
        return str(self.get_gains())
    
class DeterministicVariable(Variable):

    """
    Initialize a Deterministic Variable.

    Args:
        variable_name (str): The name of the variable.
        variables_left_gains (dict): Dictionary of variable names to their left gains.
        variables_right_gains (dict): Dictionary of variable names to their right gains.
    """

    def __init__(self, variable_name: str, variables_left_gains: Dict[str, float], variables_right_gains: Dict[str, float]):
        self.left_gain = variables_left_gains[variable_name]
        self.right_gain = variables_right_gains[variable_name]

# Parameters
gap: float = 20000 # dual gap
num_repeats: int = 100 # number of expriments for each dual gaap
total_sb_calls: int = 0 
total_tree_size: int = 0
total_more: int = 0 

# Path to the text file to read variables and their dual gains
file_path = './trento1output.txt'

# Function to read variables and their gains 
def read_variables_and_gains(file_path: str): #-> Tuple[List[str], Dict[str, float], Dict[str, float]]:
   
    """
    Read variables and their gains from a text file.

    Args:
        file_path (str): The path to the text file.

    Returns:
        variables (list): List of variable names suitable for branching at the root node.
        variables_left_gains (dict): Dictionary of variable names and their left gains.
        variables_right_gains (dict): Dictionary of variable names and their right gains.
    """

    variables_left_gains = {}
    variables_right_gains = {}
    variables_goemean_gains = {}
    variables = []

    try:
        # Define markers for the target line and the specific line
        target_marker = "[branch_fullstrong.c:558] debug: Execlp method of fullstrong branching"
        #specific_line_marker = re.compile(r'\s*\d+\.\d+s\|\s+1\s+\|\s+2\s+\|\s+')
        specific_line_marker = re.compile(r'\s*\d+s\|\s+1\s+\|\s+2\s+\|\s+.*')

        # Read the text file
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        # Find the index of the specific line
        
        specific_line_index = None
        for index, line in enumerate(lines):
            
            if specific_line_marker.search(line):
                specific_line_index = index
                break

        # Find the last occurrence of the target line before the specific line
        if specific_line_index is not None:
            last_occurrence_index = None
            for index in range(specific_line_index - 1, -1, -1):
                if target_marker in lines[index]:
                    last_occurrence_index = index
                    

                    break

            # Extract the last occurrence of the target line
            if last_occurrence_index is not None:
                desired_block = lines[last_occurrence_index + 1:specific_line_index]
                # Extract a list of variables and their up and down dual gains
                for line in desired_block:
                    
                    if line.strip().startswith("[branch_fullstrong.c:512]"):
                        var_start_index = line.find("<") + 1
                        var_end_index = line.find(">", var_start_index)
                        var_character = line[var_start_index:var_end_index]

                        variables.append(var_character)

                        downgain_start_index = line.find("downgain=") + 9
                        downgain_end_index = line.find(",", downgain_start_index)
                        downgain_number = float(line[downgain_start_index:downgain_end_index])

                        upgain_start_index = line.find("upgain=") + 7
                        upgain_end_index = line.find(",", upgain_start_index)
                        upgain_number = float(line[upgain_start_index:upgain_end_index])

                        variables_left_gains[var_character] = upgain_number
                        variables_right_gains[var_character] = downgain_number
                        variables_goemean_gains[var_character] = math.sqrt ((upgain_number + 0.01) * (downgain_number + 0.01)) - 0.01


                        
                        
        # valid_variables = {var: gain for var, gain in variables_goemean_gains.items() if gain >= 0.0001}
        # invalid_variables = {var: gain for var, gain in variables_goemean_gains.items() if gain < 0.0001}

        # Count the number of variables with geomean gain less than 0.01
        # num_invalid_variables = len(invalid_variables)
        # percent_invalid_variables = (num_invalid_variables / len(variables_goemean_gains)) 



        return variables, variables_left_gains, variables_right_gains

    except Exception as e:
        raise Exception(f"Error reading variables and gains from the file: {str(e)}")
    

def equation(x, r, l):
    try:
        result = math.pow(x, r) - math.pow(x, r - l) - 1
        if not math.isfinite(result):  # Check if the result is finite
            raise OverflowError("Result is not finite")
        return result
    except OverflowError as e:
        # print(f"Error in equation: {e}")
        return None  # Return None if the result is out of range

def find_x(r, l, tolerance=1e-6, max_iterations=10000):
    if l > r:
        l, r = r, l
    if l < 0.00001 or r < 0.00001:
        return None
    # x_low = 2**(1/r)
    # x_high = 2**(1/l)

    try:
        x_low = 2 ** (1/r)
        x_high = 2 ** (1/l)
    except OverflowError as e:
        print(f"Error in find_x: {e}")
        return None

    for _ in range(max_iterations):
        x_mid = (x_low + x_high) / 2
        result_mid = equation(x_mid, r, l)

        if result_mid is None:
            # Skip this value if it's out of range
            continue

        if abs(x_high - x_low) < tolerance:
            return x_mid

        if result_mid > 0:
            x_high = x_mid
        else:
            x_low = x_mid

    # print("l=", l)
    # print("r=", r)

    # print("Bisection method did not converge within the maximum number of iterations.")
    return None  # Return None if the method did not converge



def expected_tree_size_phi (phi:float, param:float, zero_prob, sb_count: int) -> float:
        """
        Estimates the size of two trees:
            - The tree that uses the current best variable.
            - The tree that uses the best variable given one additional strong branching iteration.
        To estimate both of these tree sizes, we use phi^gap.
        We suppose that the phi of the next variable follows a distribution that is mixed:
            - With probability zero_prob it is very large.
            - With probability (1-zero_prob) it is reasonable.

        @param phi: the growth rate of the tree.
        @param param: the parameter of the Pareto distribution.
        @param zero_prob: the probability of a dual gain close to 0, and therefore of a phi to be too large.
                          this informs 
        @param sb_count: the number of strong branching iterations so far.
        @return the estimates of the current and expected tree sizes.
        """

        if 0.1 * phi >=1:
            min_m = 0.8 * phi

        else:
            min_m = 1

        current_treesize_no_SB = phi ** gap
        current_treesize_with_SB = current_treesize_no_SB + 2 * sb_count 

        # we use the integral for phi=1 to the current phi to compute the average estimated treesize
        # with a new variable using a pareto distribution with parameter param.
        # this gives us the treesize if phi is reasonable.
        next_treesize_no_SB =  ((((param * ((min_m) ** param))/ (gap - param)) * (phi ** (gap - param))) - (((param * ((min_m) ** param))/ (gap - param)) * (min_m ** (gap- param)) )) 

        # the probability that next_treesize_no_SB is smaller than current_treesize_no_SB
        # if phi is reasonable.
        cdf = 1 - (1 / phi) ** param

        next_treesize_with_SB = zero_prob * (current_treesize_no_SB)\
             + (1- zero_prob) * (cdf * next_treesize_no_SB + (1 - cdf) * current_treesize_no_SB)\
                 + 2 * (sb_count + 1)

        # If max_variable_gain exceeds or equals the proposed gap, 
        # no further improvement in tree size is expected.
        return current_treesize_with_SB, next_treesize_with_SB

# Function for strong branching
def strong_branching(variables: List[str]) -> Tuple[int, float, int, str, float, float]:
    
    """
    Perform strong branching to select the best variable to branch on.

    Args:
        variables (list): List of variable names to choose from.

    Returns:
        sb_count (int): The number of strong branching iterations performed.
        max_variable_gain (float): The maximum geometric mean gain among the selected variables.
        more (int): A flag indicating whether more strong branching than the initial sample was performed.
        best_variable (str): The name of the best variable to branch on.
        upgain (float): The gain of the best variable when branching up.
        downgain (float): The gain of the best variable when branching down.

    Note:
        Initially, this function obtains a sample of variables, finding their dual gains, 
        and computes the current tree size based on the best variable found so far, and
        it estimates the distribution parameters. 
        It then estimates the expected tree size if one more iteration of strong branching is performed. 
        Strong branching continues until the expected tree size becomes smaller than the current tree size.
        The criteria for stopping strong branching is determined based on these expected and current tree sizes.
    """

    max_variable_gain = -math.inf
    min_variable_gain = +math.inf
    min_phi = +math.inf
    best_variable = None 
    current_treesize = 1 # tree size for the current SB iteration
    next_treesize = 0 # estimated tree size for the next SB iteration
    sb_count = 0
    changed_variable_count = 0
    more = 0
    queue = []
    phi_values = []
    nzerogains = 0
    
    # initial sample to estimate distribution parameters for subsequent decisions.
    # the strategy employed to obtain this initial sample is similar to the SCIP strategy 
    # for deciding when to stop strong branching.
    while (changed_variable_count < 12 and variables) or  best_variable is None or len(phi_values) < 12:
        
        # Randomly select a variable from the availible options.
        chosen_variable = random.choice (variables)
        # Calculate left and right dual gains for the chosen variable.
        random_variable_gains = DeterministicVariable(chosen_variable, variables_left_gains = variables_left_gains, variables_right_gains = variables_right_gains)
        # Increase the strong branching iteration count.
        sb_count += 1
        variables.remove(chosen_variable) # Remove the chosen variable from available options.


        # if random_variable_gains.left_gain == 0:
        #     random_variable_gains.left_gain = 0.000001

        # if random_variable_gains.right_gain == 0:
        #     random_variable_gains.right_gain = 0.000001
        
        # To deal with each variable easier we calculate the geometric mean of left and right dual gains.
        random_variable_geomean = math.sqrt ((random_variable_gains.left_gain + 0.01) * (random_variable_gains.right_gain + 0.01)) - 0.01

        phi = find_x(random_variable_gains.right_gain, random_variable_gains.left_gain)

        if phi is None:
            nzerogains +=1
            changed_variable_count += 1


        if phi is not None:

            if random_variable_geomean < 0.00001:
                nzerogains += 1

            else:
                phi_values. append (phi)


            if phi < min_phi:
                min_phi = phi
                best_variable = chosen_variable
                changed_variable_count = 0

            elif phi == min_phi:
                changed_variable_count += 0.5
    
            
            else:
                changed_variable_count += 1


    ln_sum = sum(math.log(phi) for phi in phi_values)
    param = len(phi_values)/ (ln_sum - len(phi_values) * math.log(min_phi))
    zero_prob = nzerogains / (nzerogains + len (phi_values))


    current_treesize, next_treesize =  expected_tree_size_phi (min_phi, param,zero_prob, sb_count)


    # compute the expected size of the tree with one more strong branching

    # Based on the expected tree size and current tree size, determine when to stop strong branching.
    while (next_treesize < current_treesize and variables):
        
        #print("Start of Strong Branching")
        #ngirihgn

        
        more = 1
        sb_count += 1 
  
        chosen_variable = random.choice(variables)  
        random_variable_gains = DeterministicVariable(chosen_variable, variables_left_gains = variables_left_gains, variables_right_gains = variables_right_gains)

        variables.remove(chosen_variable)


        #assert random_variable.left_gain >= lower_bound
        #assert random_variable.left_gain <= upper_bound
        # if random_variable_gains.left_gain == 0:
        #     random_variable_gains.left_gain = 0.000001

        # if random_variable_gains.right_gain == 0:
        #     random_variable_gains.right_gain = 0.000001

        
        random_variable_geomean = math.sqrt ((random_variable_gains.left_gain + 0.01) * (random_variable_gains.right_gain + 0.01)) - 0.01
        phi = find_x(random_variable_gains.right_gain, random_variable_gains.left_gain)
        
        if phi is None:
            nzerogains +=1

        if phi is not None:


            if random_variable_geomean < 0.00001:
                nzerogains +=1
            
            else:
                phi_values. append (phi)
        
            if phi < min_phi:
                min_phi = phi
                best_variable = chosen_variable
                
        ln_sum = sum(math.log(phi) for phi in phi_values)
        param = len(phi_values)/ (ln_sum - len(phi_values) * math.log(min_phi))
        zero_prob = nzerogains / (nzerogains + len (phi_values))

        current_treesize, next_treesize =  expected_tree_size_phi (min_phi, param, zero_prob, sb_count)
        # size = compute_treesize (random_variable_gains.right_gain, random_variable_gains.left_gain) 



        #assert upper_bound >= max_variable_gain, "The upper bound of variables is {} but we found a variable with gain {}".format(upper_bound, max_variable_gain)
 
    # find left and right dual gains for the best variable after strong branching.
    right_left_gainvariable = DeterministicVariable(best_variable, variables_left_gains = variables_left_gains, variables_right_gains = variables_right_gains)

                 
    return sb_count, max_variable_gain, more, best_variable,  right_left_gainvariable.left_gain,  right_left_gainvariable.right_gain
    
# Function to compute the tree size based on Pierre's paper
def compute_treesize(right_gain, left_gain):
    
    """
    Compute the size of the tree based on Pierre's paper.

    Args:
        right_gain (float): Bigger gain when branching.
        left_gain (float): Smaller gain when branching.

    Returns:
        treesize (int): Size of the tree.
    """
    a = 1 #Lower bound of the range for branching on the right dual gain.
    b = math.ceil(gap/right_gain) #Upper bound of the range for branching on the right dual gain
    total_sum = 0
    for i in range(a, b + 1):
        n = i + math.ceil((gap-(i-1)* right_gain)/left_gain) - 1
        if a <= i <= b:
            total_sum += math.comb(n, i)

    treesize = 1 + 2 * total_sum  

    return treesize

# Loop for experiments
for _ in range(num_repeats):

    variables, variables_left_gains, variables_right_gains = read_variables_and_gains(file_path)
    # variables,_, _ , _= read_variables_and_gains(file_path)


    sb_calls, max_gain, more_sb, best_var, upgain, downgain = strong_branching(variables)

    # fix the bigger gain as a right gain and smaller one as a left gain 
    if upgain >= downgain:
        right_gain = upgain
        left_gain = downgain
    
    if downgain > upgain:
        right_gain = downgain
        left_gain = upgain

    total_more += more_sb
    total_sb_calls += sb_calls
    
    result = compute_treesize (right_gain, left_gain)
    tree_size = compute_treesize (right_gain, left_gain) + 2 * sb_calls
    print(" tree seperate", tree_size)
    total_tree_size += tree_size

print("Number of instances that do more strong branching after getting sample=", total_more)
print("Average of SB calls:", total_sb_calls/num_repeats )
print("Average Tree Size:", total_tree_size/ num_repeats)
print("gap", gap)
