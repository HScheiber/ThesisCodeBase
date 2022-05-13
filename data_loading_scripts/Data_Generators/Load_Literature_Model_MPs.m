function Lit = Load_Literature_Model_MPs
Lit = struct();

%% Dataset 1: Aragones and Vega - https://doi.org/10.1063/1.4745205

% TF Melting points in K
Lit.Aragones.TF.LiF.Rocksalt.mp = 1010;
Lit.Aragones.TF.LiCl.Rocksalt.mp = 780;
Lit.Aragones.TF.NaF.Rocksalt.mp = 611;
Lit.Aragones.TF.NaCl.Rocksalt.mp = 1082;
Lit.Aragones.TF.NaBr.Rocksalt.mp = 1022;
Lit.Aragones.TF.KF.Rocksalt.mp = 859;
Lit.Aragones.TF.KCl.Rocksalt.mp = 1039;
Lit.Aragones.TF.KBr.Rocksalt.mp = 1043;
Lit.Aragones.TF.RbF.Rocksalt.mp = 996;
Lit.Aragones.TF.RbCl.Rocksalt.mp = 1092;
Lit.Aragones.TF.RbBr.Rocksalt.mp = 1047;

% TF Errors in melting points [K]
Lit.Aragones.TF.LiF.Rocksalt.dmp = 15;
Lit.Aragones.TF.LiCl.Rocksalt.dmp = 15;
Lit.Aragones.TF.NaF.Rocksalt.dmp = 15;
Lit.Aragones.TF.NaCl.Rocksalt.dmp = 13;
Lit.Aragones.TF.NaBr.Rocksalt.dmp = 15;
Lit.Aragones.TF.KF.Rocksalt.dmp = 15;
Lit.Aragones.TF.KCl.Rocksalt.dmp = 15;
Lit.Aragones.TF.KBr.Rocksalt.dmp = 15;
Lit.Aragones.TF.RbF.Rocksalt.dmp = 15;
Lit.Aragones.TF.RbCl.Rocksalt.dmp = 15;
Lit.Aragones.TF.RbBr.Rocksalt.dmp = 15;

% JC-SPC/E Melting points in K
Lit.Aragones.JC.NaCl.Rocksalt.mp = 1285;

% JC-SPC/E Errors in melting points [K]
Lit.Aragones.JC.NaCl.Rocksalt.dmp = 5;

% Smith-Dang Melting points in K
Lit.Aragones.JCSD.NaCl.Rocksalt.mp = 1327;

% Smith-Dang Errors in melting points [K]
Lit.Aragones.JCSD.NaCl.Rocksalt.dmp = 5;

%% Dataset 2: Walz and van der Spoel - https://doi.org/10.1039/C9CC06177K

% JC-TIP3P Melting points in K
Lit.Walz.JC3P.LiF.Rocksalt.mp = 1630;
Lit.Walz.JC3P.LiCl.Rocksalt.mp = 990;
Lit.Walz.JC3P.LiBr.Rocksalt.mp = 1060;
Lit.Walz.JC3P.LiI.Rocksalt.mp = 959;
Lit.Walz.JC3P.NaF.Rocksalt.mp = 1516;
Lit.Walz.JC3P.NaCl.Rocksalt.mp = 1324;
Lit.Walz.JC3P.NaBr.Rocksalt.mp = 1270;
Lit.Walz.JC3P.NaI.Rocksalt.mp = 1134;
Lit.Walz.JC3P.KF.Rocksalt.mp = 1368;
Lit.Walz.JC3P.KCl.Rocksalt.mp = 1230;
Lit.Walz.JC3P.KBr.Rocksalt.mp = 1207;
Lit.Walz.JC3P.KI.Rocksalt.mp = 1129;
Lit.Walz.JC3P.RbF.Rocksalt.mp = 1343;
Lit.Walz.JC3P.RbCl.Rocksalt.mp = 1191;
Lit.Walz.JC3P.RbBr.Rocksalt.mp = 1178;
Lit.Walz.JC3P.RbI.Rocksalt.mp = 1107;
Lit.Walz.JC3P.CsF.Rocksalt.mp = 1330;
Lit.Walz.JC3P.CsCl.Rocksalt.mp = 1155;

% JC-TIP3P Errors in melting points [K]
Lit.Walz.JC3P.LiF.Rocksalt.dmp = 25;
Lit.Walz.JC3P.LiCl.Rocksalt.dmp = 25;
Lit.Walz.JC3P.LiBr.Rocksalt.dmp = 25;
Lit.Walz.JC3P.LiI.Rocksalt.dmp = 25;
Lit.Walz.JC3P.NaF.Rocksalt.dmp = 25;
Lit.Walz.JC3P.NaCl.Rocksalt.dmp = 25;
Lit.Walz.JC3P.NaBr.Rocksalt.dmp = 25;
Lit.Walz.JC3P.NaI.Rocksalt.dmp = 25;
Lit.Walz.JC3P.KF.Rocksalt.dmp = 25;
Lit.Walz.JC3P.KCl.Rocksalt.dmp = 25;
Lit.Walz.JC3P.KBr.Rocksalt.dmp = 25;
Lit.Walz.JC3P.KI.Rocksalt.dmp = 25;
Lit.Walz.JC3P.RbF.Rocksalt.dmp = 25;
Lit.Walz.JC3P.RbCl.Rocksalt.dmp = 25;
Lit.Walz.JC3P.RbBr.Rocksalt.dmp = 25;
Lit.Walz.JC3P.RbI.Rocksalt.dmp = 25;
Lit.Walz.JC3P.CsF.Rocksalt.dmp = 25;
Lit.Walz.JC3P.CsCl.Rocksalt.dmp = 25;

% WBK Melting points in K
Lit.Walz.WBK.LiF.Rocksalt.mp = 1155;
Lit.Walz.WBK.LiCl.Rocksalt.mp = 1053;
Lit.Walz.WBK.LiBr.Rocksalt.mp = 998;
Lit.Walz.WBK.LiI.Rocksalt.mp = 909;
Lit.Walz.WBK.NaF.Rocksalt.mp = 1178;
Lit.Walz.WBK.NaCl.Rocksalt.mp = 1074;
Lit.Walz.WBK.NaBr.Rocksalt.mp = 1045;
Lit.Walz.WBK.NaI.Rocksalt.mp = 971;
Lit.Walz.WBK.KF.Rocksalt.mp = 1081;
Lit.Walz.WBK.KCl.Rocksalt.mp = 1018;
Lit.Walz.WBK.KBr.Rocksalt.mp = 994;
Lit.Walz.WBK.KI.Rocksalt.mp = 929;
Lit.Walz.WBK.RbF.Rocksalt.mp = 1018;
Lit.Walz.WBK.RbCl.Rocksalt.mp = 991;
Lit.Walz.WBK.RbBr.Rocksalt.mp = 953;
Lit.Walz.WBK.RbI.Rocksalt.mp = 907;
Lit.Walz.WBK.CsF.Rocksalt.mp = 1005;
Lit.Walz.WBK.CsCl.Rocksalt.mp = 968;

% WBK Errors in melting points [K]
Lit.Walz.WBK.LiF.Rocksalt.dmp = 25;
Lit.Walz.WBK.LiCl.Rocksalt.dmp = 25;
Lit.Walz.WBK.LiBr.Rocksalt.dmp = 25;
Lit.Walz.WBK.LiI.Rocksalt.dmp = 25;
Lit.Walz.WBK.NaF.Rocksalt.dmp = 25;
Lit.Walz.WBK.NaCl.Rocksalt.dmp = 25;
Lit.Walz.WBK.NaBr.Rocksalt.dmp = 25;
Lit.Walz.WBK.NaI.Rocksalt.dmp = 25;
Lit.Walz.WBK.KF.Rocksalt.dmp = 25;
Lit.Walz.WBK.KCl.Rocksalt.dmp = 25;
Lit.Walz.WBK.KBr.Rocksalt.dmp = 25;
Lit.Walz.WBK.KI.Rocksalt.dmp = 25;
Lit.Walz.WBK.RbF.Rocksalt.dmp = 25;
Lit.Walz.WBK.RbCl.Rocksalt.dmp = 25;
Lit.Walz.WBK.RbBr.Rocksalt.dmp = 25;
Lit.Walz.WBK.RbI.Rocksalt.dmp = 25;
Lit.Walz.WBK.CsF.Rocksalt.dmp = 25;
Lit.Walz.WBK.CsCl.Rocksalt.dmp = 25;

%% Dataset 3: DeFever et al. - https://doi.org/10.1063/5.0012253
% TF Melting points in K
Lit.DeFever.TF.LiCl.Rocksalt.mp = 777;
Lit.DeFever.TF.NaCl.Rocksalt.mp = 1081;
Lit.DeFever.TF.KCl.Rocksalt.mp = 1038;
Lit.DeFever.TF.RbCl.Rocksalt.mp = 1091;

% TF Errors in melting points [K]
Lit.DeFever.TF.LiCl.Rocksalt.dmp = 3;
Lit.DeFever.TF.NaCl.Rocksalt.dmp = 2;
Lit.DeFever.TF.KCl.Rocksalt.dmp = 2;
Lit.DeFever.TF.RbCl.Rocksalt.dmp = 2;

% PIM Melting points in K
Lit.DeFever.PIM.LiCl.Rocksalt.mp = 1025;
Lit.DeFever.PIM.NaCl.Rocksalt.mp = 1140;
Lit.DeFever.PIM.KCl.Rocksalt.mp = 1053;
Lit.DeFever.PIM.RbCl.Rocksalt.mp = 1035;

% PIM Errors in melting points [K]
Lit.DeFever.PIM.LiCl.Rocksalt.dmp = 5;
Lit.DeFever.PIM.NaCl.Rocksalt.dmp = 10;
Lit.DeFever.PIM.KCl.Rocksalt.dmp = 3;
Lit.DeFever.PIM.RbCl.Rocksalt.dmp = 10;
end