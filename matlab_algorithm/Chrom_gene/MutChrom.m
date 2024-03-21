function chrom_new = MutChrom(chrom, mut, N, N_chrom, chrom_range, t, iter)
for i = 1:N %%N�Ǹ���������Ҳ����ÿһ���ж���ͷ����
    for j = 1:N_chrom  %N_chrom��Ⱦɫ��ڵ����������м���Ⱦɫ��
        mut_rand = rand; %�������һ������������Ȼ��Ļ���ͻ�䣬Ȼ���ø�ֵ�������Ƿ����ͻ�䡣
        if mut_rand <=mut  %mut����ͻ����ʣ�������ͻ�����ֵ�����С��0.2�Ļ���ͻ�������ֵ�Ž��л���ͻ�䴦�����߲�����ͻ�䴦��
            mut_pm = rand; %���ӻ��Ǽ���
            mut_num = rand*(1-t/iter)^2;
            if mut_pm<=0.5
                chrom(i, j)= chrom(i, j)*(1-mut_num);
            else
                chrom(i, j)= chrom(i, j)*(1+mut_num);
            end
            chrom(i, j) = IfOut(chrom(i, j), chrom_range(:, j)); %�����Ƿ�Խ��
        end
    end
end
chrom_new = chrom;%%�ѱ��촦�����Ľ�������¾�����

