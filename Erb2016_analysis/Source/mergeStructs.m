s1.a = '1';
s1.b = '2';
s1.c = '5';

s2.c = '3';
s2.d = '4';

s1orig = s1;
s1a = s1;
s1b = s1;
s1c = s1;

% overwrites duplicate fields in s1
f = fieldnames(s2);
for i = 1:length(f)
  s1a = setfield(s1, f{i}, getfield(s2, f{i}));
end

% same but with dynamic field names - works in R14 and up I believe.
f = fieldnames(s2);
for i = 1:length(f)
    s1b.(f{i}) = s2.(f{i});
end

% somewhat cryptic code - also overwrites duplicate fields in s1
s1 = rmfield(s1,intersect(fieldnames(s1),fieldnames(s2)));
v1 = [fieldnames(s1) struct2cell(s1)];
v2 = [fieldnames(s2) struct2cell(s2)];
v3 = [reshape(v1.',1,[]) reshape(v2.',1,[])];
s1c = struct(v3{:});
clear v1 v2 v3

%restore original s1 for clarity
s1 = s1orig;

%clear intermediate variables
clear s1orig f i

display(s1a)
display(s1b)
display(s1c)