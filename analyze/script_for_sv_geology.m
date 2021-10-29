for portion_str={"0.10", "0.20", "0.30", "0.40", "0.50"};
    for repeat_str={"1", "2", "3"};
        for algoname={"midas", "ns/global_deg_1.0000", "es/global_deg_min_0.0000", "ns/global_deg_0.0000", "tihs", "rw/rw_c_1", "ff/ff_c_0.51_0.20", "mgs/add_degree", "mgs/exchange_degree", "mgs/remove_degree", "mgs/add_avg", "mgs/exchange_avg", "mgs/remove_avg", "midas_grid_ablation", "maxdegree", "avgdegree", "midas_ns"};
            fname = algoname{1} + '/coauth-MAG-Geology-full_' + portion_str{1} + '/' + repeat_str{1}
            [fileID_input, msg] = fopen('../Hypergraph_Sampling_cpp/results/'+ algoname{1} + '/coauth-MAG-Geology-full_' + portion_str{1} + '/' + repeat_str{1} + '/sampled_graph.txt', 'r');
            [fileID, msg] = fopen('input/' + algoname{1} + '/coauth-MAG-Geology-full_' + portion_str{1} + '/' + repeat_str{1} + '/icd_row.txt','r');
            if (fileID_input < 0) || (fileID < 0);
                disp(msg);
                continue;
            end;
            fileID2 = fopen('input/' + algoname{1} + '/coauth-MAG-Geology-full_' + portion_str{1} + '/' + repeat_str{1} + '/icd_col.txt','r');
            fileID3 = fopen('input/'+ algoname{1}+ '/coauth-MAG-Geology-full_'+ portion_str{1}+ '/' + repeat_str{1} + '/dim.txt','r');
            formatSpec = '%f';
            R = fscanf(fileID, formatSpec);                       
            C = fscanf(fileID2, formatSpec);
            Dim = fscanf(fileID3, formatSpec);
            S = sparse(R,C,1,Dim(1), Dim(2)); 
            sv = svds(S, 300);
            sz = size(sv,1);
            dlmwrite('output/'+ algoname{1} + '/coauth-MAG-Geology-full_' + portion_str{1} + '/' + repeat_str{1} + '/singular_values_full.txt',sv(:),'newline','pc')
        end;
    end;
end;