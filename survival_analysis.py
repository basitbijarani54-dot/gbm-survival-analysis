import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

def run_merged_survival_analysis(clinical_path, mutation_path):
    print("1. Parsing custom cBioPortal survival data...")
    # Custom parser to bypass the ParserError in the text file
    with open(clinical_path, 'r') as f:
        lines = f.readlines()
        
    data = []
    for line in lines:
        parts = [p.strip() for p in line.strip('\n').split('\t')]
        # Identify valid patient data rows
        if len(parts) >= 6 and parts[0] != 'Case ID' and parts[0] != '':
            data.append({
                'Case ID': parts[0],
                'Status': parts[3],
                'OS_MONTHS': parts[5]
            })
            
    df_clinical = pd.DataFrame(data)
    if df_clinical.empty:
        raise ValueError("No data could be extracted from the survival file.")
    print(f"Loaded {len(df_clinical)} clinical records.")
    
    print("2. Loading mutation data...")
    df_mut = pd.read_csv(mutation_path, sep='\t')
    mut_id_col = next((c for c in df_mut.columns if 'SAMPLE' in c.upper() or 'PATIENT' in c.upper()), None)
    
    # Ensure IDs are strings
    mutated_ids = [str(mid) for mid in df_mut[mut_id_col].unique()] if mut_id_col else []
    print(f"Found {len(mutated_ids)} patients with IDH1 mutations.")
    
    print("3. Merging clinical and mutation data using substring matching...")
    # Tag patients dynamically using substring matching to handle TCGA ID discrepancies
    def tag_patient(case_id):
        for mid in mutated_ids:
            if str(case_id) in mid or mid in str(case_id):
                return 'IDH1 Mutated'
        return 'IDH1 Wild-Type'

    df_clinical['Group'] = df_clinical['Case ID'].apply(tag_patient)
    
    # Clean Survival Data
    df_clinical['OS_MONTHS'] = pd.to_numeric(df_clinical['OS_MONTHS'], errors='coerce')
    df_clinical = df_clinical.dropna(subset=['OS_MONTHS', 'Status'])
    
    # Encode Status (1 = Deceased, 0 = Living/Censored)
    df_clinical['Event'] = df_clinical['Status'].astype(str).apply(
        lambda x: 1 if any(term in x.upper() for term in ['1', 'DECEASED', 'DEAD', 'TRUE']) else 0
    )
    
    print("4. Generating Kaplan-Meier plot...")
    kmf = KaplanMeierFitter()
    plt.figure(figsize=(10, 6))
    
    groups = ['IDH1 Mutated', 'IDH1 Wild-Type']
    colors = {'IDH1 Mutated': 'darkred', 'IDH1 Wild-Type': 'navy'}
    
    for g in groups:
        sub_df = df_clinical[df_clinical['Group'] == g]
        if not sub_df.empty:
            kmf.fit(durations=sub_df['OS_MONTHS'], event_observed=sub_df['Event'], label=f"{g} (n={len(sub_df)})")
            kmf.plot_survival_function(color=colors[g], linewidth=2)
        
    # Log-rank Statistical Test
    g1 = df_clinical[df_clinical['Group'] == 'IDH1 Mutated']
    g2 = df_clinical[df_clinical['Group'] == 'IDH1 Wild-Type']
    
    if not g1.empty and not g2.empty:
        p_val = logrank_test(g1['OS_MONTHS'], g2['OS_MONTHS'], 
                             event_observed_A=g1['Event'], event_observed_B=g2['Event']).p_value
        p_text = f'Log-rank p-value: {p_val:.4e}'
    else:
        p_text = ''
        
    # Format the Graph
    plt.title('Kaplan-Meier Survival Analysis: IDH1 in Glioblastoma (Merged Data)', fontsize=13, fontweight='bold')
    plt.xlabel('Overall Survival (Months)', fontsize=11)
    plt.ylabel('Survival Probability', fontsize=11)
    
    if p_text:
        plt.text(0.05, 0.10, p_text, transform=plt.gca().transAxes, 
                 fontsize=11, bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray'))
                 
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    plt.savefig('survival_curve_merged.png', dpi=300)
    print("\nSUCCESS: Mutation and Clinical data successfully merged!")
    print("Analysis complete. Saved to 'survival_curve_merged.png'.")
    plt.show()

if __name__ == '__main__':
    run_merged_survival_analysis('survival_data.txt', 'idh1_mutations.tsv')