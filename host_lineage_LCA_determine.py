import pandas as pd
import sys

'''example
argv[1]
UViGs	Phylum	Class	Order	Family	Genus	Species
Virus_001	P_1	C_1	O_1	F_1	G_1	S_1
Virus_001	P_1	C_1	O_1	F_1	G_2	S_2
Virus_001	P_1	C_1	O_1	F_1	G_2	S_3
Virus_002	P_1	C_2	O_2	F_2	G_3	S_4
Virus_002	P_1	C_2	O_2	F_3	G_4	S_5
Virus_003	P_2	C_3	O_3	F_4	G_5	S_6
Virus_003	P_2	C_3	O_4	F_5	G_6	S_7
...

argv[2]
UViGs	Phylum	Class	Order	Family	Genus	Species
Virus_001	P_1	C_1	O_1	F_1	-	-
Virus_002	P_1	C_2	O_2	-	-	-
Virus_003	P_2	C_3	-	-	-	-

'''

def find_last_common_lineage(input_file, output_file):

    df = pd.read_csv(input_file, sep='\t')
    
    columns = ['UViGs', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    result_list = []
    
    # Cycle through each UViG
    for uvi_gs, group in df.groupby('UViGs'):
        # initialize a dictionary to store the final lineage
        lineage = {col: None for col in columns}
        lineage['UViGs'] = uvi_gs
        
        # initialize LCA
        common_lineage = None
        
        # Cycle through the taxonomic ranks of each potential host
        for _, row in group.iterrows():
			
            current_lineage = row.dropna().tolist()[1:]
            
            if common_lineage is None:
                common_lineage = current_lineage
            else:
                # Update the LC-lineage to the LCA of the current row
                common_lineage = [common_lineage[i] if i < len(common_lineage) and i < len(current_lineage) and common_lineage[i] == current_lineage[i] else None for i in range(len(columns) - 1)]
                
                # If the LC-lineage no longer has any common nodes, stop the comparison
                if not any(common_lineage):
                    break
        
        # Update dictionary
        for i, col in enumerate(columns[1:], 1): 
            if common_lineage and i <= len(common_lineage):
                lineage[col] = common_lineage[i-1]
        
        result_list.append(lineage)
    
    result_df = pd.concat([pd.DataFrame([lineage]) for lineage in result_list], ignore_index=True)
    
    result_df = result_df.fillna('-')
    
    result_df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    print("Copyright: Kaiyang Zheng, OUC, zhengkaiyang@stu.ouc.edu.cn")
    find_last_common_lineage(input_file, output_file)
