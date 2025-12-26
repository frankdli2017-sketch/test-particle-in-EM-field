%% set f in multi points 

   close all; clear; clc; 
    
%% add BN
% 3D streamline in B field

%     close all;clear;clc
    % 12-22 22:29 CS event
    B0 = 0; % nT
    B01 = -13; B02 = -13; 
    BG = 13.54;
    BH0 = 10;
    Bn = 3.0;
    N1 = -15; N2 = 15; N3 = -13; N4 = 13; L3 = N3*2; L4 = N4*2; % km
    a1 = 15; a2 = 15; a3 = 15; a4 = 15; a5 = a3*2; a6 = a4*2; % km
    De = 2.6; % km
        
    tic
    % figure
    figure(20)
    xSize=850; ySize=600;
    set(gcf,'Position',[100 100 xSize ySize]);
    set(gca,'Position',[0.10,0.14,0.78,0.78]);
    
%     XLim = [-60, 60];
%     YLim = [-30, 30];
%     ZLim = [-40, 80];
%     XLim = [-60, 60]*2;
%     YLim = [-30, 30]*2;
%     ZLim = [-40, 80]*2;
    XLim = [-60, 60];
    YLim = [-30, 30];
    ZLim = [-20, 80];

    % quiver3
    Lmax0 = 320*2; Mmax0 = 240*2; Nmax0 = 160*2; 
    L2_v0 = linspace(-Lmax0,Lmax0,81);
    Y2_v0 = linspace(-Nmax0,Nmax0,81);
%     L2_v0 = -Lmax0:16:Lmax0; % km
%     Y2_v0 = -Nmax0:8:Nmax0; % km, N
%     L2_v0 = -160:16:160; % km
%     Y2_v0 = -80:8:80; % km, N
%     Z2_v0 = -120:24:120; % km, M
    Z2_v0 = 0*De; % km, M
    
    [LL2_v0,YY2_v0,ZZ2_v0] = meshgrid(L2_v0,Y2_v0,Z2_v0);
    BL22_v0 = B0 + B01*tanh((YY2_v0-N1)/a1) + B02*tanh((YY2_v0-N2)/a2);
    BY22_v0 = YY2_v0*0 + Bn;
    BG22_v0 = YY2_v0*0 + BG;
    BH22_v0 = BH0*(exp(-(YY2_v0-N3).^2/a3^2)-exp(-(YY2_v0-N4).^2/a4^2)).*(exp(-(LL2_v0-L3).^2/a5^2)-exp(-(LL2_v0-L4).^2/a6^2));
    BZ22_v0 = BG22_v0 + BH22_v0;
    De = 2.6; % km
    LL2_v0 = LL2_v0/De; ZZ2_v0 = ZZ2_v0/De; YY2_v0 = YY2_v0/De;
    q = quiver3(LL2_v0, YY2_v0,ZZ2_v0, BL22_v0,BY22_v0,BZ22_v0,  'Color', [0.2 0.2 0.2]);
    q.LineWidth = 1.2;
    q.MaxHeadSize = 0.2;
    
%     [startx,starty,startz] = meshgrid([10],[-20:5:20],[0]);
%     s = streamline(LL2_v0, YY2_v0,ZZ2_v0, BL22_v0, BY22_v0,BZ22_v0,  startx,starty,startz);
    
    % stream line 
%     Lmax0 = 320; Mmax0 = 240; Nmax0 = 160; 
    Lmax0 = 320*2; Mmax0 = 240*2; Nmax0 = 160*2; 
    L2_v1 = linspace(-Lmax0,Lmax0,201);
    Y2_v1 = linspace(-Nmax0,Nmax0,201);
    Z2_v1 = linspace(-Mmax0,Mmax0,201);
    
%     L2_v1 = -160:0.8:160; % km
%     Z2_v1 = -120:0.6:120; % km
%     Y2_v1 = -80:0.4:80; % km
%     L2_v1 = -Lmax0:0.8:Lmax0; % km
%     Y2_v1 = -Nmax0:0.4:Nmax0; % km
%     Z2_v1 = -Mmax0:0.6:Mmax0; % km
    
    [LL2_v1,YY2_v1,ZZ2_v1] = meshgrid(L2_v1,Y2_v1,Z2_v1);
    
    BL22_v1 = B0 + B01*tanh((YY2_v1-N1)/a1) + B02*tanh((YY2_v1-N2)/a2);
    BY22_v1 = YY2_v1*0 + Bn;
    BG22_v1 = YY2_v1*0 + BG;
    BH22_v1 = BH0*(exp(-(YY2_v1-N3).^2/a3^2)-exp(-(YY2_v1-N4).^2/a4^2)).*(exp(-(LL2_v1-L3).^2/a5^2)-exp(-(LL2_v1-L4).^2/a6^2));
    BZ22_v1 = BG22_v1 + BH22_v1;
    De = 2.6; % km
    LL2_v1 = LL2_v1/De; YY2_v1 = YY2_v1/De; ZZ2_v1 = ZZ2_v1/De;
    
%     Ltmp1 = -100:10:50; Ntmp1 = -30:5:-5;
    Ltmp1 = -150:10:50; Ntmp1 = -15; Mtmp1 = -30;
    [startx,starty,startz] = meshgrid(Ltmp1,Ntmp1,Mtmp1);
    s2 = streamline(LL2_v1, YY2_v1, ZZ2_v1, BL22_v1, BY22_v1, BZ22_v1, startx,starty,startz);
    for i = 1:length(s2)
        s2(i).Color = [0.8 0 0];
        s2(i).LineStyle = '-';
        s2(i).LineWidth = 1.5;
    end
    
%     Ltmp1 = -100:10:100; Ntmp1 = 0;
%     [startx,starty,startz] = meshgrid(Ltmp1,Ntmp1,[0]);
%     s = streamline(LL2_v1, YY2_v1, ZZ2_v1, BL22_v1, BY22_v1, BZ22_v1, startx,starty,startz);
%     for i = 1:length(s)
%         s(i).Color = [0 0.6 0];
%         s(i).LineStyle = '-';
%         s(i).LineWidth = 1.5;
%     end
%     
% %     Ltmp1 = -50:10:100; Ntmp1 = 5:5:30;
%     Ltmp1 = -50:10:100; Ntmp1 = 4;
%     [startx,starty,startz] = meshgrid(Ltmp1,Ntmp1,[0]);
%     s = streamline(LL2_v1, YY2_v1, ZZ2_v1, BL22_v1, BY22_v1, BZ22_v1, startx,starty,startz);
%     for i = 1:length(s)
%         s(i).Color = [0 0 0.8];
%         s(i).LineStyle = '-';
%         s(i).LineWidth = 1.5;
%     end
    
    toc
    
    xlabel(gca,'L (d_e)');
    ylabel(gca,'N (d_e)');
    zlabel(gca,'M (d_e)');
%     title(gca,['R_{0,LMN} = [',num2str(P1(1)),' ,',num2str(P1(2)),' ,',num2str(P1(3)),' ] de,    V_{0,LMN} = [',num2str(P1(4)),' ,',num2str(P1(5)),' ,',num2str(P1(6)),'] 10^3 km/s'],'fontsize',16)
    axis(gca, 'equal');
    %
    set(gca,'xlim',XLim); 
    set(gca,'ylim',YLim);
    set(gca,'Zlim',ZLim);
    set(gca,'YDir','reverse')
    set(gcf,'color','w')
    set(gca,'fontsize',18);
%     view(gca, [-15 30])
    view(3)

    
%% multi particles
%     close all
    tic
    % 12-22 22:29 CS event
    B0 = 0; % nT
    B01 = -13; B02 = -13; 
    BG = 13.54;
    BH0 = 10;
    Bn = 3.0;
    N1 = -15; N2 = 15; N3 = -13; N4 = 13; L3 = N3*2; L4 = N4*2; % km
    a1 = 15; a2 = 15; a3 = 15; a4 = 15; a5 = a3*2; a6 = a4*2; % km
    
    N = -80:0.4:80; % km
    L = -160:0.8:160; % km
    [LL,NN] = meshgrid(L,N);
    
    BL2 = B0 + B01*tanh((NN-N1)/a1) + B02*tanh((NN-N2)/a2);
    BN2 = NN*0 + Bn;
    BG2 = NN*0 + BG;
    BH2 = BH0*(exp(-(NN-N3).^2/a3^2)-exp(-(NN-N4).^2/a4^2)).*(exp(-(LL-L3).^2/a5^2)-exp(-(LL-L4).^2/a6^2));
    BM2 = BG2 + BH2;
    Bt2 = sqrt(BL2.^2+BM2.^2+BN2.^2);
    De = 2.6; % km
    LL = LL/De; NN = NN/De;
    Xl0 = 140/De; Zl0 = 70/De;
    
    % ----------- test particle ----------------
    % 5.1. parameters
    Dem = 2600;                    % electron inertial length; [m]
    mp = 0;                         % electron;
    % Axis-X:L, Axis-Y:-N, Axis-Z:M, 
    B0 = 0e-9; % nT
    B01 = -13e-9; B02 = -13e-9; % two step Harris current sheet B0
    BG = 13.54e-9; % Guide field 
    BH0 = 10e-9;
    Bn = 3.0e-9;
    N1 = -15; N2 = 15; N3 = -13; N4 = 13; L3 = N3*2; L4 = N4*2; % km
    DNa1 = 15; DNa2 = 15; DNa3 = 15; DNa4 = 15; DNa5 = DNa3*2; DNa6 = DNa4*2; % km, Harris current sheet half thickness 
    NZ = [N1,N2,N3,N4,L3,L4]*1e3; % m
    DNZ = [DNa1,DNa2,DNa3,DNa4,DNa5,DNa6]*1e3; % m
    
    EX0 = 0.;
    EY0 = 0.; % 3.0e-3;                      % [V/m]
    EZ0 = 0.;
    
    dsweep = 100;                    % sweep length: 10 m;
    trange = [0,0.080]; % s
    nsweep = 8000;
    % polyfit parameters [BL -- EN];
    pp = [0    0    0];     
    
    % figure
    figure(21)
    xSize=850; ySize=600;
    set(gcf,'Position',[100 100 xSize ySize]);
    set(gca,'Position',[0.10,0.14,0.78,0.78]);
    
%     XLim = [-Xl0, Xl0];
%     YLim = [-Zl0, Zl0];
%     ZLim = [-30,30]; 
%     XLim = [-60, 60];
%     YLim = [-30, 30];
%     ZLim = [-40, 40];

    XLim = [-60, 60];
    YLim = [-30, 30];
    ZLim = [-20, 80];
    
    color_bg_up = [220 220 220]/256;
    color_bg_down = [230 245 242]/256;
    color_curve = [194, 80, 124]/256;
    color_curve_pro = [0, 0, 0,100]/256;
    
    % --------------- set parameters ----------------
    Vtt = 10; % 1e6 m/s
%     ePA = 10:20:170; % degree

%     ePA  = 180;
%     Vtt = 5:5:15; % 1e6 m/s

%     ePA = 0:10:90; % degree
%     Ltmp = -10; Mtmp = 10; Ntmp = -8; % de

    ePA = 0:10:80; % degree
    ptmp = [-60.2 -14.8 -24.7]; % xyz

%     ePA = 90:10:180; % degree
% %     ptmp = [-60.0 14.9 109.4]; % xyz
%     ptmp = [-70.0 14.9 109.4]; % xyz

        
    Ltmp = ptmp(1); Mtmp = ptmp(3); Ntmp = ptmp(2); % de

%     ePA = 90:10:180; % degree
%     ePA = 90:10:270; % degree
%     Ltmp = 7; Mtmp = 70; Ntmp = 6; % de
    
    
    
    % get color map
    clrs = get_clrs(length(ePA));
%     clrs = get_clrs(length(Vtt));

%     Vtt = 5:5:40; % 1e6 m/s
%     ePA = 25; % degree
%     % get color map
%     clrs = get_clrs(length(Vtt));
    
    hold on
    hd1 = surf(LL,NN,NN*0+0,Bt2);
    shading('flat');
    colormap(gca,othercolor('BuDRd_18')); 
    h = colorbar;
    h.Label.String = 'B_t (nT)';
    h.Position = [0.8725 0.3100 0.0174 0.4150];
    
    
    % test particle parameter

%     Ltmp = 10; Mtmp = 40; Ntmp = -4; % de
%     Ltmp = 30; Mtmp = 70; Ntmp = -4; % de
%     Ltmp = 40; Mtmp = 60; Ntmp = -4; % de
%     Ltmp = 30; Mtmp = 60; Ntmp = -4; % de
%     Ltmp = 30; Mtmp = 70; Ntmp = -4; % de
%     Ltmp = -10; Mtmp = 10; Ntmp = -4; % de
    
    Bvectmp = get_Bvec(Ltmp,Mtmp,Ntmp,Dem,B0,B01,B02,BG,BH0, Bn, NZ, DNZ);
    % VLtmp =-20; VMtmp = 0; VNtmp =0;
    tmpdir = [0,0,1];
%     tmpdir = [0,1,0];
    tmpdir2 = cross(tmpdir,Bvectmp);
    tmpdir2 = tmpdir2/norm(tmpdir2);
    
    %
    PAtot = [];
    tttot = [];
    RRxtot = [];
    RRytot = [];
    RRztot = [];
    VVxtot = [];
    VVytot = [];
    VVztot = [];
    
%     for kk = 1:length(Vtt)
% 
%         tmpdir3 = Bvectmp*cosd(ePA)+tmpdir2*sind(ePA);
%         VLtmp = Vtt(kk)*tmpdir3(1);
%         VMtmp = Vtt(kk)*tmpdir3(3);
%         VNtmp = Vtt(kk)*tmpdir3(2);
        
    for kk = 1:length(ePA)

        tmpdir3 = Bvectmp*cosd(ePA(kk))+tmpdir2*sind(ePA(kk));
        VLtmp = Vtt*tmpdir3(1);
        VMtmp = Vtt*tmpdir3(3);
        VNtmp = Vtt*tmpdir3(2);

        P1 = [Ltmp,Mtmp,Ntmp, VLtmp,VMtmp,VNtmp];
        [PA1, tt1, RR1, VV1] = get_testparticle(P1,mp,Dem,B0,B01,B02,BG,BH0, Bn, EX0, EY0, EZ0, NZ, DNZ, pp, dsweep,trange,nsweep);
        plot3(RR1(1, 1)/Dem, RR1(1, 2)/Dem, RR1(1, 3)/Dem, 'o','MarkerSize', 9,'LineWidth', 1.8, 'color','k');  
        plot3(RR1(:, 1)/Dem, RR1(:, 2)/Dem, RR1(:, 3)/Dem, 'LineWidth', 2.0, 'color', clrs(kk,:)); 
%         plot3(RR1(:, 1)/Dem, RR1(:, 2)/Dem, RR1(:, 3)/Dem, 'LineWidth', 2.0, 'color', [0.1+abs(VLtmp)/50,0,0.9-abs(VLtmp)/50]); 
        PAtot = [PAtot,real(PA1)];
        tttot = [tttot,real(tt1)];
        RRxtot = [RRxtot,RR1(:,1)];
        RRytot = [RRytot,RR1(:,2)];
        RRztot = [RRztot,RR1(:,3)];
        VVxtot = [VVxtot,VV1(:,1)];
        VVytot = [VVytot,VV1(:,2)];
        VVztot = [VVztot,VV1(:,3)];
    
    end
    
    % quiver3
    q = quiver3(LL2_v0, YY2_v0,ZZ2_v0*0+0, BL22_v0,BY22_v0,BZ22_v0,  'Color', [0.2 0.2 0.2]);
    q.LineWidth = 1.2;
    q.MaxHeadSize = 0.2;
    
    % stream line
%     Ltmp1 = -100:10:50; Ntmp1 = -10; Mtmp1 = 0;
    Ltmp1 = -150:10:50; Ntmp1 = -16; Mtmp1 = -30;
    [startx,starty,startz] = meshgrid(Ltmp1,Ntmp1,Mtmp1);
    s1 = streamline(LL2_v1, YY2_v1, ZZ2_v1, BL22_v1, BY22_v1, BZ22_v1, startx,starty,startz);
    for i = 1:length(s1)
%         s(i).Color = [0.8 0 0];
        s1(i).Color = 'm';
        s1(i).LineStyle = '-';
        s1(i).LineWidth = 1.8;
    end
    
    % stream line
    Ltmp1 = -100:10:50; Ntmp1 = -4; Mtmp1 = 0;
    [startx,starty,startz] = meshgrid(Ltmp1,Ntmp1,Mtmp1);
    s2 = streamline(LL2_v1, YY2_v1, ZZ2_v1, BL22_v1, BY22_v1, BZ22_v1, startx,starty,startz);
    for i = 1:length(s2)
        s2(i).Color = [0 0 0];
        s2(i).LineStyle = '-';
        s2(i).LineWidth = 1.8;
    end
    
%     % stream line
%     Ltmp1 = -100:10:100; Ntmp1 = 0; Mtmp1 = 0;
%     [startx,starty,startz] = meshgrid(Ltmp1,Ntmp1,Mtmp1);
%     s = streamline(LL2_v1, YY2_v1, ZZ2_v1, BL22_v1, BY22_v1, BZ22_v1, startx,starty,startz);
%     for i = 1:length(s)
%         s(i).Color = [0 0.4 0];
%         s(i).LineStyle = '-';
%         s(i).LineWidth = 1.8;
%     end
    
%     % stream line
% %     Ltmp1 = -50:10:100; Ntmp1 = 5:5:30;
%     Ltmp1 = -50:10:100; Ntmp1 = 10;  Mtmp1 = 0;
%     [startx,starty,startz] = meshgrid(Ltmp1,Ntmp1,Mtmp1 );
%     s = streamline(LL2_v1, YY2_v1, ZZ2_v1, BL22_v1, BY22_v1, BZ22_v1, startx,starty,startz);
%     for i = 1:length(s)
%         s(i).Color = [0 0 0.8];
%         s(i).LineStyle = '-';
%         s(i).LineWidth = 1.8;
%     end
    
    
    % background color
    fill3([XLim(1) XLim(2) XLim(2) XLim(1)], [YLim(1), YLim(1), YLim(1), YLim(1)], ...
        [ZLim(1), ZLim(1), ZLim(2), ZLim(2)], color_bg_up,'FaceAlpha',0.4);       % left up
    fill3([XLim(2) XLim(2) XLim(2) XLim(2)], [YLim(2), YLim(1), YLim(1), YLim(2)], ...
        [ZLim(1), ZLim(1), ZLim(2), ZLim(2)], color_bg_up,'FaceAlpha',0.4);       % right up
    fill3([XLim(1) XLim(2) XLim(2) XLim(1)], [YLim(1), YLim(1), YLim(2), YLim(2)], ...
        [ZLim(1), ZLim(1), ZLim(1), ZLim(1)], color_bg_up,'FaceAlpha',0.4);       % down
    
%     plot3([XLim(1) XLim(2)],[-Ntmp,-Ntmp],[ZLim(1), ZLim(1)],'color','b','Linestyle','--','Linewidth',2.5)
%     plot3([XLim(1) XLim(2)],[-10,-10],[ZLim(1), ZLim(1)],'color','m','Linestyle','--','Linewidth',2.5)
    
    hold off
    xlabel(gca,'L (d_e)');
    ylabel(gca,'N (d_e)');
    zlabel(gca,'M (d_e)');
%     title(gca,['R_{0,LMN} = [',num2str(P1(1)),' ,',num2str(P1(2)),' ,',num2str(P1(3)),' ] de,    V_{0,LMN} = [',num2str(P1(4)),' ,',num2str(P1(5)),' ,',num2str(P1(6)),'] 10^3 km/s'],'fontsize',16)
    axis(gca, 'equal');
    %
    set(gca,'xlim',XLim); 
    set(gca,'ylim',YLim);
    set(gca,'Zlim',ZLim);
    grid on
    caxis(gca, [7.5,30])
    set(gca,'YDir','reverse')
    set(gcf,'color','w')
    set(gca,'fontsize',18);
    view(gca, [-15 30])
    toc
    

    %% pitch angle

    figure(31)
    xSize=850; ySize=250;
    set(gcf,'Position',[100 400 xSize ySize]);
    set(gca,'Position',[0.10,0.22,0.84,0.73]);
    
    tic
    
    hold on
    for kk = 1:length(ePA)
        
%         tmpdir3 = Bvectmp*cosd(ePA(kk))+tmpdir2*sind(ePA(kk));
%         VLtmp = Vtt*tmpdir3(1);
%         VMtmp = Vtt*tmpdir3(3);
%         VNtmp = Vtt*tmpdir3(2);
%         
% %         VLtmp = Vtt*cosd(ePA(kk))*Bvectmp(1);
% %         VMtmp = Vtt*cosd(ePA(kk))*Bvectmp(3);
% %         VNtmp = Vtt*sind(ePA(kk));
% 
%         P1 = [Ltmp,Mtmp,Ntmp, VLtmp,VMtmp,VNtmp];
%         [PA1, tt1, RR1, VV1] = get_testparticle(P1,mp,Dem,B0,B01,B02,BG,BH0, Bn, EX0, EY0, EZ0, NZ, DNZ, pp, dsweep,trange,nsweep);
%         
        plot(tt1*1e3, PAtot(:,kk), 'color',clrs(kk,:),'linestyle','-','linewidth',1.2)

    end
    
    plot([0,100], [90,90], 'color','k','linestyle','--','linewidth',1.2)

    hold off
    toc
    
    xlabel('T (ms)');
    ylabel('PA (º)');
    set(gca,'xlim',[0,80]);
    set(gca,'ylim',[0 180]);
    grid on
    set(gca,'fontsize',17)
    set(gca,'box','on')
    set(gcf,'color','w')
    hd1 = gca;
    hd1.LineWidth = 1.0;
    
    %% particle trajectory
    
%     figure(41)
%     xSize=850; ySize=250;
%     set(gcf,'Position',[100 400 xSize ySize]);
%     set(gca,'Position',[0.10,0.22,0.84,0.73]);
    
    
    npanel1 = 3; 
    h = irf_plot(npanel1,'newfigure'); 
    xSize=850; ySize=750;
    set(gcf,'Position',[100 100 xSize ySize]);
    
    set(h(1),'Position',[0.12 0.12 0.20 0.60]);
    set(h(2),'Position',[0.32 0.12 0.60 0.60]);
    set(h(3),'Position',[0.32 0.72 0.60 0.20]);
    
    
    ip = 1;
    hold(h(ip),'on')
    for ii = 6:15 % 1:length(s1)
        plot(h(ip),s1(ii).YData, s1(ii).XData, 'color','m','linestyle','-','linewidth',1.0)
    end
%     for ii = 1:length(s2)
%         plot(h(ip),s2(ii).YData, s2(ii).XData, 'color','k','linestyle','-','linewidth',1.2)
%     end
    for kk = 1:length(ePA)
        plot(h(ip),RRytot(:,kk)/Dem, RRxtot(:,kk)/Dem, 'color',clrs(kk,:),'linestyle','-','linewidth',1.5)
    end
        plot(h(ip),RRytot(1,kk)/Dem, RRxtot(1,kk)/Dem, 'o','MarkerSize', 10,'LineWidth', 1.8, 'color','k');  
        
    hold(h(ip),'off')
    xlabel(h(ip),'N (de)');
    ylabel(h(ip),'L (de)');
    set(h(ip),'xlim',[-32,32]);
    set(h(ip),'ylim',[-80 80]);
    
    ip = 2;
    hold(h(ip),'on')
 
    for ii = 6:15 % 1:length(s1)
        plot(h(ip),s1(ii).ZData, s1(ii).XData, 'color','m','linestyle','-','linewidth',1.0)
    end
    for ii = 9 % 1:length(s1)
        plot(h(ip),s1(ii).ZData, s1(ii).XData, 'color','m','linestyle','-','linewidth',2.5)
    end
%     for ii = 1:length(s2)
%         plot(h(ip),s2(ii).ZData, s2(ii).XData, 'color','k','linestyle','-','linewidth',1.2)
%     end
    for kk = 1:length(ePA)
        plot(h(ip),RRztot(:,kk)/Dem, RRxtot(:,kk)/Dem, 'color',clrs(kk,:),'linestyle','-','linewidth',1.5)
    end
        plot(h(ip),RRztot(1,kk)/Dem, RRxtot(1,kk)/Dem, 'o','MarkerSize', 10,'LineWidth', 1.8, 'color','k'); 
        
    hold(h(ip),'off')
    xlabel(h(ip),'M (de)');
    ylabel(h(ip),' ');
    set(h(ip),'xlim',[-40,120]);
    set(h(ip),'ylim',[-80 80]);
    set(h(ip),'yticklabel',' ');
    
    ip = 3;
    hold(h(ip),'on')
    for kk = 1:length(ePA)
        plot(h(ip),RRztot(:,kk)/Dem, RRytot(:,kk)/Dem, 'color',clrs(kk,:),'linestyle','-','linewidth',1.2)
    end
        plot(h(ip),RRztot(1,kk)/Dem, RRytot(1,kk)/Dem, 'o','MarkerSize', 10,'LineWidth', 1.8, 'color','k');  
%     for ii = 1:length(s1)
%         plot(h(ip),s1(ii).ZData, s1(ii).YData, 'color','m','linestyle','-','linewidth',1.2)
%     end
        
    hold(h(ip),'off')
    xlabel(h(ip),' ');
    ylabel(h(ip),'N (de)');
    set(h(ip),'xlim',[-40,120]);
    set(h(ip),'ylim',[-32 32]);
    set(h(ip),'xticklabel',' ');
    
    for ii = 1:3
        grid(h(ii),'on')
        set(h(ii),'fontsize',17)
        set(h(ii),'box','on')
        h(ii).LineWidth = 1.0;
    end
    
    %
%     pos0 = [s1(9).XData;s1(9).YData;s1(9).ZData]';
    
    
    

    %% PA ~ N
    
    figure(51)
    xSize=850; ySize=250;
    set(gcf,'Position',[100 400 xSize ySize]);
    set(gca,'Position',[0.10,0.22,0.84,0.73]);
    
    tic
    hold on
    for kk = 1:9 %[1,3,5,7]
        plot(RRytot(:,kk)/Dem, PAtot(:,kk), 'color',clrs(kk,:),'linestyle','-','linewidth',1.2)
    end
    plot([-100,100], [90,90], 'color','k','linestyle','--','linewidth',1.2)
    hold off
    toc

    xlabel('N (de)');
    ylabel('PA (º)');
    set(gca,'xlim',[-20,20]);
    set(gca,'ylim',[0 180]);
    grid on
    set(gca,'fontsize',17)
    set(gca,'box','on')
    set(gcf,'color','w')
    hd1 = gca;
    hd1.LineWidth = 1.0;
    
    %%
    
    
    
    
    %% -------------------------------------------
    function result = get_clrs(ll)
        clrs0 = colormap(jet);
        s1=21; s2=241;
        
        if ll == 1
            result = clrs0(s1,:);
        else
            ds = floor((s2-s1)/(ll-1));
            ss = s1:ds:s2;
            result = clrs0(ss,:);
        end
        
    end
    
    
    
    
    function Bvec = get_Bvec(L,M,N,Dem,B0,B01,B02,BG,BH0, Bn, NZ, DNZ)
        
        BL2 = B0 + B01*tanh((N*Dem-NZ(1))/DNZ(1)) + B02*tanh((N*Dem-NZ(2))/DNZ(2));
        BH2 = BH0*(exp(-(N*Dem-NZ(3)).^2/DNZ(3)^2)-exp(-(N*Dem-NZ(4)).^2/DNZ(4)^2)).*(exp(-(L*Dem-NZ(5)).^2/DNZ(5)^2)-exp(-(L*Dem-NZ(6)).^2/DNZ(6)^2));
        BM2 = BG + BH2;
        
        Bsweep = [BL2, Bn, BM2];
        Bvec = Bsweep./norm(Bsweep);
        
    end
    
    %%
    function [PA1, tt1, RR1, VV1] = get_testparticle(P1,mp,Dem,B0,B01,B02,BG,BH0, Bn, EX0, EY0, EZ0, NZ, DNZ, pp, dsweep,trange,nsweep)
        R0 = [P1(1), P1(3), P1(2)]*Dem;               
        V0 = [P1(4), P1(6), P1(5)]*1e6;  % 1e6 m/s --> m/s  

%         trange1 = [0 0.010];                 % [s];
        trange1 = trange;                 % [s];
        Tarrow1 = 1;                        % time forward [1: forward; -1 backward];   
        % test particle
        [tt1, RR1, VV1] = leapfrog_lorentz_solver04(mp, R0, V0, trange1, Tarrow1, B0,B01,B02,BG,BH0, Bn, EX0, EY0, EZ0, NZ, DNZ, pp, dsweep,nsweep); 
        [Atheta1,Aphi1,PA1] = pitchangle(RR1, VV1, B0,B01,B02,BG,BH0, Bn, NZ, DNZ);
        
    end
    

    function results = get_testparticle_v2(P1,mp,Dem,B0,B01,B02,BG,BH0, Bn, EX0, EY0, EZ0, NZ, DNZ, pp, dsweep,trange,nsweep)
    
        R0 = [P1(1), P1(3), P1(2)]*Dem;               
        V0 = [P1(4), P1(6), P1(5)]*1e6; % 1e6 m/s --> m/s  

%         trange1 = [0 0.010];                 % [s];
        trange1 = trange;                 % [s];
        Tarrow1 = 1;                        % time forward [1: forward; -1 backward];   
        % test particle
        [tt1, RR1, VV1] = leapfrog_lorentz_solver04(mp, R0, V0, trange1, Tarrow1, B0,B01,B02,BG,BH0, Bn, EX0, EY0, EZ0, NZ, DNZ, pp, dsweep,nsweep); 
        [Atheta1,Aphi1,PA1] = pitchangle(RR1, VV1, B0,B01,B02,BG,BH0, Bn, NZ, DNZ);
        
        results = cell(6,1);
        results{1} = Atheta1; results{2} = Aphi1; results{3} = PA1; results{4} = tt1; results{5} = RR1; results{6} = VV1;
        
    end

%% ---------------------------- Leapfrog Solver ---------------------------
    function [tt, RR, VV] = leapfrog_lorentz_solver04(mp, R0, V0, trange, Tarrow, B0,B01,B02,BG,BH0, Bn, EX0, EY0, EZ0, NZ, DNZ, pp, dsweep,nsweep)
%     tic;
    units = irf_units; 
    if mp == 0 
        Mp = units.me;
        Qp = units.e * (-1);
    else
        Mp = units.mp * mp;
        Qp = units.e * mp;
    end
    Vminus = V0;
    RR = R0;
    VV = V0;
    % make Time reversal parameter;
    if Tarrow > 0.5
        Trev = 1;
        Vrev = 1;
        Brev = 1;
        Erev = 1;
        disp(['*NOTE*: Time forward ', num2str(trange(1)), ' --------> ', num2str(trange(2))]);
    elseif Tarrow < -0.5
        Trev = -1;
        Vrev = -1;
        Brev = -1;
        Erev = 1;   
        disp(['*NOTE*: ', num2str(trange(2)), ' <-------- ', num2str(trange(1)), ' Time backward']);
    end    
    
%  load E & B field
    r0 = R0;
 
    % T-Symmetry for initia velocity & time
    Vminus = Vminus * Vrev;
    trange = trange * Trev;
    t0 = trange(1);
    tt = trange(1);

    while (t0 <=trange(2))
        % 1. compute EB field

        BL2 = B0 + B01*tanh((r0(2)-NZ(1))/DNZ(1)) + B02*tanh((r0(2)-NZ(2))/DNZ(2));
        BH2 = BH0*(exp(-(r0(2)-NZ(3)).^2/DNZ(3)^2)-exp(-(r0(2)-NZ(4)).^2/DNZ(4)^2)).*(exp(-(r0(1)-NZ(5)).^2/DNZ(5)^2)-exp(-(r0(1)-NZ(6)).^2/DNZ(6)^2));
        BM2 = BG + BH2;
        
        Bsweep = [BL2, Bn, BM2];
        
        Esweep = [EX0, EY0, EZ0];
        Esweep(3) = polyval(pp, Bsweep(1)*1e9) * 1e-3;
        
        % DT: delte Time
            % fce = abs(Qp) * BX0/Mp/2/pi;
            % DT = dt * 1/fce;
%             Vamp_tmp = norm(Vminus);
%             DT = dsweep / Vamp_tmp;
%             DT=0.03/3e3;
            DT = trange(2)/nsweep;
        % for T-Symmetry
            Bsweep = Brev * Bsweep;
            Esweep = Erev * Esweep;
    
        % Leapfrog scheam
        vminus = Vminus + DT * Qp * Esweep / 2 / Mp;        % Eqn. (8-5-6)
        tvec = DT/2 * Qp * Bsweep / Mp;                     % Eqn. (8-5-12)
        svec = DT * Qp * norm(Bsweep) / Mp /(1+(norm(tvec))^2) * Bsweep/norm(Bsweep); % Eqn. (8-5-13)
        avec = vminus + cross(vminus, tvec);                % Eqn. (8-5-14)
        vplus = vminus + cross(avec, svec);                 % Eqn. (8-5-15)
        Vplus = vplus + DT * Qp/2/Mp * Esweep;              % Eqn. (8-5-5)
        
        % one step forward
        tt = [tt; tt(end)+DT];
        rr = r0 + DT * Vplus;
        RR = [RR; rr];
        VV = [VV; Vplus];
        
        % prepare for next step
        t0 = tt(end);
        r0 = rr;
        Vminus = Vplus;
    end
    VV = VV * Tarrow;               % if Tarrow = -1; output VV is along time forward direction.
%     toc;
    end    
    
%% ------------------------------ Pitch Angle -----------------------------
    function [Atheta,Aphi,PA] = pitchangle(RR, VV, B0,B01,B02,BG,BH0, Bn, NZ, DNZ)
    
        % 1. basic
        ntt = length(RR);
        PA = zeros(ntt, 1);
        Atheta = zeros(ntt, 1);
        Aphi = zeros(ntt, 1);
        for ii = 1 : ntt
            rr = RR(ii, :);
            vv = VV(ii, :);

            BL2 = B0 + B01*tanh((rr(2)-NZ(1))/DNZ(1)) + B02*tanh((rr(2)-NZ(2))/DNZ(2));
            BH2 = BH0*(exp(-(rr(2)-NZ(3)).^2/DNZ(3)^2)-exp(-(rr(2)-NZ(4)).^2/DNZ(4)^2)).*(exp(-(rr(1)-NZ(5)).^2/DNZ(5)^2)-exp(-(rr(1)-NZ(6)).^2/DNZ(6)^2));
            BM2 = BG + BH2;

            Btmp = [BL2, Bn, BM2];

            PA(ii) = acosd(dot(vv, Btmp)/norm(vv)/norm(Btmp));
            Atheta(ii) = acosd(vv(3)./norm(vv)); % tht_angle: polar angle
            if vv(2)>0
                 Aphi(ii) = acosd(vv(1)./sqrt(vv(1)^2+vv(2)^2));
            else
                 Aphi(ii) = 360-acosd(vv(1)./sqrt(vv(1)^2+vv(2)^2));
            end
            
        end
    end
    
    %%
    
    
    
    
    
    
