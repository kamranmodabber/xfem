function [B,J0] = xfemBmatrix1(pt,elemType,e,enrich_node,...
                              Qq,noke_tarak,zaviye,noke_gere)
global node element

sctr = element(e,:);
nn   = length(sctr);
[N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
invJ0 = inv(J0);
dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
Gpt = N' * node(sctr,:);                  % GP in global coord, used
%faghat hesab karde

% Bfem is always computed
Bfem = zeros(3,2*nn);
Bfem(1,1:2:2*nn)  = dNdx(:,1)' ;
Bfem(2,2:2:2*nn)  = dNdx(:,2)' ;
Bfem(3,1:2:2*nn)  = dNdx(:,2)' ;
Bfem(3,2:2:2*nn)  = dNdx(:,1)' ;

% Switch between non-enriched and enriched elements
if ( any(enrich_node(sctr)) == 0 ) % Non-enriched elements
    B = Bfem;
else                               % Enriched elements
    Bxfem = [] ;
    % loop on nodes, check node is enriched ...
    for in = 1 : nn

        
        if ( enrich_node(sctr(in)) == 1)     % H(x) enriched node
      
            % Enrichment function, H(x) at global Gauss point
           
            dist = signed_distance1(Qq,Gpt);
            Hgp  = heaviside(dist);
            % Enrichment function, H(x) at node "in"
            dist = signed_distance1(Qq,node(sctr(in),:));
            Hi   = heaviside(dist);
            % Bxfem at node "in"
            BI_enr = [dNdx(in,1)*(Hgp - Hi) 0 ;
                0 dNdx(in,2)*(Hgp - Hi) ;
                dNdx(in,2)*(Hgp - Hi) dNdx(in,1)*(Hgp - Hi)];
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
        elseif ( enrich_node(sctr(in)) == 2) % B(x) enriched node
 if noke_gere(sctr(in),2)==1
            
              xTip=noke_tarak(1,2:3);
  
              alpha=zaviye(1,1);
              QT  =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            else
          %xCr=Qq(end-1:end,:) ; 
          xTip=noke_tarak(2,2:3);
  
         alpha=zaviye(end,2);
         QT  =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
            end
            
            % compute branch functions at Gauss point
            xp    = QT*(Gpt-xTip)';           % local coordinates
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            if ( theta > pi | theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            [Br,dBdx,dBdy] = branch(r,theta,alpha);

            % compute branch functions at node "in"
            xp    = QT*(node(sctr(in),:)-xTip)';
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
            theta = atan2(xp(2),xp(1));
            if ( theta > pi | theta < -pi)
                disp (['something wrong with angle ',num2str(thet)]);
            end
            [BrI] = branch_node(r,theta);

            % composants of Benr matrix

            aa = dNdx(in,1)*(Br(1)-BrI(1)) + N(in)*dBdx(1) ;
            bb = dNdx(in,2)*(Br(1)-BrI(1)) + N(in)*dBdy(1) ;
            B1_enr = [aa 0 ; 0 bb ; bb aa];

            aa = dNdx(in,1)*(Br(2)-BrI(2)) + N(in)*dBdx(2) ;
            bb = dNdx(in,2)*(Br(2)-BrI(2)) + N(in)*dBdy(2) ;
            B2_enr = [aa 0 ; 0 bb ; bb aa];

            aa = dNdx(in,1)*(Br(3)-BrI(3)) + N(in)*dBdx(3) ;
            bb = dNdx(in,2)*(Br(3)-BrI(3)) + N(in)*dBdy(3) ;
            B3_enr = [aa 0 ; 0 bb ; bb aa];

            aa = dNdx(in,1)*(Br(4)-BrI(4)) + N(in)*dBdx(4) ;
            bb = dNdx(in,2)*(Br(4)-BrI(4)) + N(in)*dBdy(4) ;
            B4_enr = [aa 0 ; 0 bb ; bb aa];

            BI_enr = [B1_enr B2_enr B3_enr B4_enr];
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
        end
    end          % end of loop on nodes
    % B matrix
    B = [ Bfem Bxfem ];
    clear Bfem; clear Bxfem;
end              % end of switch between enriched and non-enriched elements

