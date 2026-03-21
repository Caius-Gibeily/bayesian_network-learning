function new_s=my_py2mat(a,b,c)
n=py.py2mat.test(a,b,c);
s=struct(n);
field_list=fieldnames(s);
rows=double(py.len(s.(field_list{1}))); %size and numel can't acquire the dimension of a python dictionary, that's why we call py.len()
new_array=cell(rows,numel(field_list));
for i=1:numel(field_list)
    temp=s.(field_list{i});
    for j=1:rows
        if isa(temp{j-1},'py.str')
            new_array{j,i}=char(temp{j-1});
        elseif isa(temp{j-1},'py.int')
            new_array{j,i}=double(temp{j-1});
        else 
            new_array{j,i}=temp{j-1};
        end
    end
     s.(field_list{i})=new_array(:,i);
end
new_s=s;

% credit to:
% Ajpaezm: https://www.mathworks.com/matlabcentral/answers/371814-transform-a-dict-or-any-python-data-structure-to-struct-or-cell-array-in-matlab