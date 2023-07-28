function  [c_Line,c_Poly,c_LWPoly,c_Cir,c_Arc,c_Poi] = read_dxf(nomArch)
%%  Read dxf file with different type like point line polyline and 3dface
% written by Prof Meng from Hohai University 
% 
fId = fopen(nomArch);
c_ValAsoc = textscan(fId,'%d%s','Delimiter','\n');
fclose(fId);
% Code Group Matrix
m_GrCode = c_ValAsoc{1};
% Associated value String Cell
c_ValAsoc = c_ValAsoc{2};
%% Entities
m_PosCero = find(m_GrCode==0);
%Is searched by (0,SECTION),(2,ENTITIES)
indInSecEnt = strmatch('ENTITIES',c_ValAsoc(m_PosCero(1:end-1)+1),'exact');
%(0,ENDSEC)
m_indFinSecEnt = strmatch('ENDSEC',c_ValAsoc(m_PosCero(indInSecEnt:end)),'exact');
% Entities Position
m_PosCero = m_PosCero(indInSecEnt:indInSecEnt-1+m_indFinSecEnt(1));
I=0;
    c_Line = cell(1,2);
    c_Poly = cell(1,2);
    c_LWPoly = cell(1,2);
    c_Cir = cell(1,2);
    c_Arc = cell(1,2);
    c_Poi = cell(1,2);
    c_Face = cell(1,2);
    %
    iLine = 1;
    iPoly = 0;
    iLWPoly = 1;
    iCir = 1;  
    iArc = 1;
    iPoi = 1;
    iFace = 1;
    % Loop on the Entities
for iEnt = 1:length(m_PosCero)-2
    m_GrCodeEnt = m_GrCode(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
    c_ValAsocEnt = c_ValAsoc(m_PosCero(iEnt+1):m_PosCero(iEnt+2)-1);
    nomEnt = c_ValAsocEnt{1};
%     disp(nomEnt);
    switch nomEnt 
            case 'LINE'
                % (Xi,Yi,Zi,Xj,Yj,Zj) start and end points
                c_Line{iLine,1} = [str2double(c_ValAsocEnt(m_GrCodeEnt==10)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==20)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==30)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==11)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==21)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==31))];
                % Layer
                c_Line(iLine,2) =c_ValAsocEnt(m_GrCodeEnt==8);
                % Color
                %if no exist is ByLayer (256)
                %c_Line{iLine,3} = str2double(f_ValGrCode(62,m_GrCodeEnt,c_ValAsocEnt));
                % XData
                %c_Line(iLine,4) = f_XData(GroupCode,'XDataName',m_GrCodeEnt,c_ValAsocEnt);
                % Add properties
                %
                iLine = iLine+1;  
        case 'POLYLINE'
            TMP=zeros(1,2);j=0;iPoly=iPoly+1;
        case 'VERTEX'
            j=j+1;
            TMP(j,:)=[str2double(c_ValAsocEnt(m_GrCodeEnt==10)),str2double(c_ValAsocEnt(m_GrCodeEnt==20))];
            %DAT{I}=TMP;
            if str2double(c_ValAsocEnt(m_GrCodeEnt==10))==1
                TMP(j+1,:)=TMP(1,:);
            end
            c_Poly{iPoly,1}=TMP;
            c_Poly(iPoly,2)=c_ValAsocEnt(m_GrCodeEnt==8);
            
        case 'LWPOLYLINE'
                % (X,Y) vertexs
                %Is not take into account the budge (group code 42, arc in the polyline).
                m_Coord = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                        str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt))];
                if strcmp(f_ValGrCode(70,m_GrCodeEnt,c_ValAsocEnt),'1')&&...
                        any(m_Coord(1,:)~=m_Coord(end,:))
                    %Close polyline
                    c_LWPoly{iLWPoly,1} = [m_Coord;m_Coord(1,:)];
                else
                    c_LWPoly{iLWPoly,1} = m_Coord;
                end
                % Layer
                c_LWPoly(iLWPoly,2) = c_ValAsocEnt(m_GrCodeEnt==8);
                % Add properties
                %
                iLWPoly = iLWPoly+1;
        case '3DFACE'
            %I=I+1;
            p1=[str2double(c_ValAsocEnt(m_GrCodeEnt==10)),str2double(c_ValAsocEnt(m_GrCodeEnt==20)),str2double(c_ValAsocEnt(m_GrCodeEnt==30))];
            p2=[str2double(c_ValAsocEnt(m_GrCodeEnt==11)),str2double(c_ValAsocEnt(m_GrCodeEnt==21)),str2double(c_ValAsocEnt(m_GrCodeEnt==31))];
            p3=[str2double(c_ValAsocEnt(m_GrCodeEnt==12)),str2double(c_ValAsocEnt(m_GrCodeEnt==22)),str2double(c_ValAsocEnt(m_GrCodeEnt==32))];
            p4=[str2double(c_ValAsocEnt(m_GrCodeEnt==13)),str2double(c_ValAsocEnt(m_GrCodeEnt==23)),str2double(c_ValAsocEnt(m_GrCodeEnt==33))];
            if sum(abs(p3-p4))<0.01
                c_Face{iFace,1}=[p1;p2;p3];
                %disp('ok');
            else
                c_Face{iFace,1}=[p1;p2;p3;p4];
            end
            c_Face(iFace,2) = c_ValAsocEnt(m_GrCodeEnt==8);
            iFace=iFace+1;
            case 'CIRCLE'
                % (X Center,Y Center,Radius)
                c_Cir{iCir,1} = [str2double(f_ValGrCode(10,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(20,m_GrCodeEnt,c_ValAsocEnt)),...
                    str2double(f_ValGrCode(40,m_GrCodeEnt,c_ValAsocEnt))];
                % Layer
                c_Cir(iCir,2) = c_ValAsocEnt(m_GrCodeEnt==8);
                % Add properties
                %
                iCir = iCir+1;
            case 'ARC'
                % (X Center,Y Center,Radius,Start angle,End angle)
                c_Arc{iArc,1} = [str2double(c_ValAsocEnt(m_GrCodeEnt==10)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==20)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==40)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==50)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==51))];
                % Layer
                c_Arc(iArc,2) = c_ValAsocEnt(m_GrCodeEnt==8);
                % Add properties
                %
                iArc = iArc+1;
            case 'POINT'
                % (X,Y,Z) Position
                c_Poi{iPoi,1} = [str2double(c_ValAsocEnt(m_GrCodeEnt==10)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==20)),...
                    str2double(c_ValAsocEnt(m_GrCodeEnt==30))];
                % Layer
                c_Poi(iPoi,2) = c_ValAsocEnt(m_GrCodeEnt==8);
                % Add properties
                %
                iPoi = iPoi+1;
            %case Add Entities
    end
end
end

