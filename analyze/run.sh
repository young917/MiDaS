dataset=("email-Enron-full" "email-Eu-full" "NDC-classes-full" "NDC-substances-full" "contact-high-school" "contact-primary-school" "tags-ask-ubuntu" "tags-math-sx" "threads-ask-ubuntu" "coauth-MAG-Geology-full" "coauth-MAG-History-full")
types=("baseline" "compete" "ablation")

for t in ${types[@]}
do
    # make csv table
    python analyze_result.py --type $t
    # analyze robustness
    python draw_figures.py --type $t --select 0
    # get figures of distributions in P1-P10
    for p in 0.1 0.2 0.3 0.4 0.5
    do
        for data in ${dataset[@]}
        do
            python draw_figures.py --type $t --select 2 --dataname $data --portion $p
        done
    done
done

# observations
python observation.py --select 0
for p in 0.10 0.20 0.30 0.40 0.50
do
    python observation.py --select 1 --portion $p
done
for data in ${dataset[@]}
do
    for p in 0.10 0.20 0.30 0.40 0.50    
    do
        python observation.py --select 2 --dataname $data --portion $p
    done
    python observation.py --select 3
    python observation.py --select 4 --dataname $data
done

# ablation study
python ablation_study.py --select 0 --algotype es --opt global_deg_min
python ablation_study.py --select 0 --algotype ns --opt global_deg
python ablation_study.py --select 0 --algotype es --opt global_deg_max
python ablation_study.py --select 0 --algotype es --opt global_deg_avg
for data in ${dataset[@]}
do
    for p in 0.10 0.20 0.30 0.40 0.50    
    do
        python ablation_study.py --select 1 --algotype ns --opt global_deg --data $data --portion $p
        python ablation_study.py --select 1 --algotype es --opt global_deg_max --data $data --portion $p
        python ablation_study.py --select 1 --algotype es --opt global_deg_avg --data $data --portion $p
        python ablation_study.py --select 2 --data $data --portion $p
    done
done
python ablation_study.py --select 3

# theoretical analysis
for data in ${dataset[@]}
do
    python theorem_plot.py --data $data
done