function plot_psd_diff(data, freq_range, lim, diffAorB, byTC, figname)

figure('Name', figname)

    switch byTC
        case 0
            if ~isempty(data) 
                switch diffAorB
                case 'A'
                    plot(freq_range, [data.psddiffA])
                case 'B' 
                    plot(freq_range, [data.psddiffB])
                end
                legend({data.code})
                title(strcat(getVarName(data), " diff ", diffAorB))
                ylim(lim)
            end


        case 1
            if ~isempty(data) 
            subplot(1,2,1)
            switch diffAorB
                case 'A'
                    plot(freq_range, [data.psddiffA_TC])
                case 'B' 
                    plot(freq_range, [data.psddiffB_TC])
            end

                legend({data.code})
                title(strcat(getVarName(data), " diff ", diffAorB, "- TC contacts"))
                ylim(lim)
            subplot(1,2,2)
            switch diffAorB
                case 'A'
                    plot(freq_range, [data.psddiffA_noTC])
                case 'B' 
                    plot(freq_range, [data.psddiffB_noTC])
            end
                legend({data.code})
                title(strcat(getVarName(data), " diff ", diffAorB, "- no TC contacts"))
                ylim(lim)
            end
    end
end
       