function out=cutstr(str,query,ev)
if nargin<3
	ev=0;
end
out='';
n1=strfind(str,[ char(39),query,char(39)]);
if ~isempty(n1) 
	
	n21=strfind(str(n1(1):length(str)),',...');
    if ~strcmp(query,'Data')
		n22=strfind(str(n1(1):length(str)),');');
	else
		n22=n21;
	end

    n2=min([n21 n22]);
	
	if ~isempty(n2) 
        out=str(n1(1)+length(query)+3:n1(1)+n2(1)-2);
		out=str(n1(1)+length(query)+3:n1(1)+n2(1)-2);
    end
	
	if strfind(out,'@')
		%out='';
	end
	
	if strfind(out,'get')
		%out;
	end	
	
	if ev==0
		if ~isempty(out)
			out=eval(out);
		end	
	end
	
end	
