
<<EOF
1 安装好指定lucaprot git目录到yaml文件中
2 本地直接读取00-results目录的l。fasta文件分析结果到l.csv目录下。
3 最终结果需要下一步整合的代码读取到同一个文件中
EOF
cd LucaProt/src/
python predict_many_samples.py \
	--fasta_file /workspace/liangpeiyu/project/hantaKG_Agent/results_cell_score/00_results/l.fasta  \
	--save_file /workspace/liangpeiyu/project/hantaKG_Agent/results_cell_score/00_results/l.fasta.csv  \
	--emb_dir ../emb/   \
	--truncation_seq_length 4096  \
	--dataset_name rdrp_40_extend  \
	--dataset_type protein     \
	--task_type binary_class     \
	--model_type sefn     \
	--time_str 20230201140320   \
	--step 100000  \
	--threshold 0.5 \
	--print_per_number 10 \
	--gpu_id -1