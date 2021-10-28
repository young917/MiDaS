dataset=("email-Eu-full" "coauth-MAG-Geology-full" "contact-primary-school")
searchspace=("0.0" "0.125" "0.1768" "0.25" "0.3536" "0.5" "0.7071" "1.0" "1.4142" "2.0" "2.8284" "4.0" "5.6569" "8.0" "11.3137" "16.0" "22.6274" "32.0" "45.2548" "64.0")

for repeat_index in 1 2 3
do
    for p in 0.1 0.2 0.3 0.4 0.5
    do
        for data in ${dataset[@]}
        do
            # MiDaS & MiDaS_Basic
            for alpha in ${searchspace[@]}
            do
                ./bin/Sampling --dataname $data --inputpath ../dataset/ --algorithm es --alpha $alpha --algo_opt global_deg_min --repeat $i --target_portion $portion
            done
            # Simple and Intuitive Approaches
            ./bin/Sampling --algorithm tihs --dataname $data --inputpath ../dataset/ --repeat $repeat_index --target_portion $portion
            ./bin/Sampling --algorithm rw --dataname $data --maxlength 1 --inputpath ../dataset/ --algo_opt rw_c --repeat $repeat_index --target_portion $portion
            ./bin/Sampling --algorithm ff --dataname $data --inputpath ../dataset/ --algo_opt ff_c --repeat $repeat_index --target_portion $portion
            ./bin/Sampling --algorithm ns --dataname $data --inputpath ../dataset/ --algo_opt global_deg --alpha 0 --repeat $repeat_index --target_portion $portion
            ./bin/Sampling --algorithm ns --dataname $data --inputpath ../dataset/ --algo_opt global_deg --alpha 1 --repeat $repeat_index --target_portion $portion
            # Metropolis Graph Sampling
            ./bin/Sampling --repeat $repeat_index --eval_opt avg --algorithm mgs --dataname $data --inputpath ../dataset/ --algo_opt remove --target_portion $portion
            ./bin/Sampling --repeat $repeat_index --eval_opt avg --algorithm mgs --dataname $data --inputpath ../dataset/ --algo_opt add --target_portion $portion
            ./bin/Sampling --repeat $repeat_index --eval_opt avg --algorithm mgs --dataname $data --inputpath ../dataset/ --algo_opt exchange --target_portion $portion
            ./bin/Sampling --repeat $repeat_index --eval_opt degree --algorithm mgs --dataname $data --inputpath ../dataset/ --algo_opt remove --target_portion $portion
            ./bin/Sampling --repeat $repeat_index --eval_opt degree --algorithm mgs --dataname $data --inputpath ../dataset/ --algo_opt add --target_portion $portion
            ./bin/Sampling --repeat $repeat_index --eval_opt degree --algorithm mgs --dataname $data --inputpath ../dataset/ --algo_opt exchange --target_portion $portion
        done
    done
done
