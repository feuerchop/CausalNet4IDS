function [data2, sym_tbl] = merge_data(data)
% data is a text-numeric mixed cell array
    data2 = [];
    sym_tbl = {};
    for i = 1:size(data, 2)
        col = data{1, i};
        if iscell(col)
            % it is a cell containing string
            symbols = {};
            v = 0;
            for j = 1:length(col)
                s = col(j);
                if isempty(find(strcmpi(s, symbols)))
                    v = v + 1;
                    symbols(length(symbols) + 1, 1) = s;
                    data2(j, i) = v;
                else
                    data2(j, i) = find(strcmpi(s, symbols));
                end
            end
            sym_tbl{size(sym_tbl, 1) + 1, 1} = i;
            sym_tbl{size(sym_tbl, 1), 2} = symbols;
        else
            data2(:, i) = col;
        end
    end
end                