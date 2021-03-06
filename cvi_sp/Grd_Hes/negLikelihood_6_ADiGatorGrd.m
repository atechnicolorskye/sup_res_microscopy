% This code was generated using ADiGator version 1.1.1
% �2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function out = negLikelihood_6_ADiGatorGrd(x,w,y,c,baseline,one,Xt,Yt)
global ADiGator_negLikelihood_6_ADiGatorGrd
if isempty(ADiGator_negLikelihood_6_ADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_negLikelihood_6_ADiGatorGrd.negLikelihood_6_ADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
cada1f1 = w(1);
cada1f2 = cada1f1*1;
cada1f3 = cada1f2/one;
cada1f4 = cada1f3*1;
cada1f5 = abs(c);
cada1f6 = cada1f5^2;
cada1f7 = Gator1Data.Data1*cada1f6;
cada1f8 = cada1f4/cada1f7;
cada1f9dx = x.dx(1);
cada1f9 = x.f(1);
cada1tempdx = cada1f9dx(Gator1Data.Index1);
cada1f10dx = -cada1tempdx;
cada1f10 = Xt - cada1f9;
cada1f11dx = 2.*cada1f10(:).^(2-1).*cada1f10dx;
cada1f11 = cada1f10.^2;
cada1f12dx = x.dx(7);
cada1f12 = x.f(7);
cada1tempdx = cada1f12dx(Gator1Data.Index2);
cada1f13dx = -cada1tempdx;
cada1f13 = Yt - cada1f12;
cada1f14dx = 2.*cada1f13(:).^(2-1).*cada1f13dx;
cada1f14 = cada1f13.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index3) = cada1f11dx;
cada1td1(Gator1Data.Index4) = cada1td1(Gator1Data.Index4) + cada1f14dx;
cada1f15dx = cada1td1;
cada1f15 = cada1f11 + cada1f14;
cada1f16dx = -cada1f15dx;
cada1f16 = uminus(cada1f15);
cada1f17 = c^2;
cada1f18 = 2*cada1f17;
cada1f19dx = cada1f16dx./cada1f18;
cada1f19 = cada1f16/cada1f18;
cada1tf1 = cada1f19(Gator1Data.Index5);
cada1f20dx = exp(cada1tf1(:)).*cada1f19dx;
cada1f20 = exp(cada1f19);
cada1f21dx = cada1f8.*cada1f20dx;
cada1f21 = cada1f8*cada1f20;
cada1f22dx = cada1f21dx;
cada1f22 = baseline + cada1f21;
cada1f23 = w(2);
cada1f24 = cada1f23*1;
cada1f25 = cada1f24/one;
cada1f26 = cada1f25*1;
cada1f27 = abs(c);
cada1f28 = cada1f27^2;
cada1f29 = Gator1Data.Data2*cada1f28;
cada1f30 = cada1f26/cada1f29;
cada1f31dx = x.dx(2);
cada1f31 = x.f(2);
cada1tempdx = cada1f31dx(Gator1Data.Index6);
cada1f32dx = -cada1tempdx;
cada1f32 = Xt - cada1f31;
cada1f33dx = 2.*cada1f32(:).^(2-1).*cada1f32dx;
cada1f33 = cada1f32.^2;
cada1f34dx = x.dx(8);
cada1f34 = x.f(8);
cada1tempdx = cada1f34dx(Gator1Data.Index7);
cada1f35dx = -cada1tempdx;
cada1f35 = Yt - cada1f34;
cada1f36dx = 2.*cada1f35(:).^(2-1).*cada1f35dx;
cada1f36 = cada1f35.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index8) = cada1f33dx;
cada1td1(Gator1Data.Index9) = cada1td1(Gator1Data.Index9) + cada1f36dx;
cada1f37dx = cada1td1;
cada1f37 = cada1f33 + cada1f36;
cada1f38dx = -cada1f37dx;
cada1f38 = uminus(cada1f37);
cada1f39 = c^2;
cada1f40 = 2*cada1f39;
cada1f41dx = cada1f38dx./cada1f40;
cada1f41 = cada1f38/cada1f40;
cada1tf1 = cada1f41(Gator1Data.Index10);
cada1f42dx = exp(cada1tf1(:)).*cada1f41dx;
cada1f42 = exp(cada1f41);
cada1f43dx = cada1f30.*cada1f42dx;
cada1f43 = cada1f30*cada1f42;
cada1td1 = zeros(4096,1);
cada1td1(Gator1Data.Index11) = cada1f22dx;
cada1td1(Gator1Data.Index12) = cada1td1(Gator1Data.Index12) + cada1f43dx;
cada1f44dx = cada1td1;
cada1f44 = cada1f22 + cada1f43;
cada1f45 = w(3);
cada1f46 = cada1f45*1;
cada1f47 = cada1f46/one;
cada1f48 = cada1f47*1;
cada1f49 = abs(c);
cada1f50 = cada1f49^2;
cada1f51 = Gator1Data.Data3*cada1f50;
cada1f52 = cada1f48/cada1f51;
cada1f53dx = x.dx(3);
cada1f53 = x.f(3);
cada1tempdx = cada1f53dx(Gator1Data.Index13);
cada1f54dx = -cada1tempdx;
cada1f54 = Xt - cada1f53;
cada1f55dx = 2.*cada1f54(:).^(2-1).*cada1f54dx;
cada1f55 = cada1f54.^2;
cada1f56dx = x.dx(9);
cada1f56 = x.f(9);
cada1tempdx = cada1f56dx(Gator1Data.Index14);
cada1f57dx = -cada1tempdx;
cada1f57 = Yt - cada1f56;
cada1f58dx = 2.*cada1f57(:).^(2-1).*cada1f57dx;
cada1f58 = cada1f57.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index15) = cada1f55dx;
cada1td1(Gator1Data.Index16) = cada1td1(Gator1Data.Index16) + cada1f58dx;
cada1f59dx = cada1td1;
cada1f59 = cada1f55 + cada1f58;
cada1f60dx = -cada1f59dx;
cada1f60 = uminus(cada1f59);
cada1f61 = c^2;
cada1f62 = 2*cada1f61;
cada1f63dx = cada1f60dx./cada1f62;
cada1f63 = cada1f60/cada1f62;
cada1tf1 = cada1f63(Gator1Data.Index17);
cada1f64dx = exp(cada1tf1(:)).*cada1f63dx;
cada1f64 = exp(cada1f63);
cada1f65dx = cada1f52.*cada1f64dx;
cada1f65 = cada1f52*cada1f64;
cada1td1 = zeros(6144,1);
cada1td1(Gator1Data.Index18) = cada1f44dx;
cada1td1(Gator1Data.Index19) = cada1td1(Gator1Data.Index19) + cada1f65dx;
cada1f66dx = cada1td1;
cada1f66 = cada1f44 + cada1f65;
cada1f67 = w(4);
cada1f68 = cada1f67*1;
cada1f69 = cada1f68/one;
cada1f70 = cada1f69*1;
cada1f71 = abs(c);
cada1f72 = cada1f71^2;
cada1f73 = Gator1Data.Data4*cada1f72;
cada1f74 = cada1f70/cada1f73;
cada1f75dx = x.dx(4);
cada1f75 = x.f(4);
cada1tempdx = cada1f75dx(Gator1Data.Index20);
cada1f76dx = -cada1tempdx;
cada1f76 = Xt - cada1f75;
cada1f77dx = 2.*cada1f76(:).^(2-1).*cada1f76dx;
cada1f77 = cada1f76.^2;
cada1f78dx = x.dx(10);
cada1f78 = x.f(10);
cada1tempdx = cada1f78dx(Gator1Data.Index21);
cada1f79dx = -cada1tempdx;
cada1f79 = Yt - cada1f78;
cada1f80dx = 2.*cada1f79(:).^(2-1).*cada1f79dx;
cada1f80 = cada1f79.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index22) = cada1f77dx;
cada1td1(Gator1Data.Index23) = cada1td1(Gator1Data.Index23) + cada1f80dx;
cada1f81dx = cada1td1;
cada1f81 = cada1f77 + cada1f80;
cada1f82dx = -cada1f81dx;
cada1f82 = uminus(cada1f81);
cada1f83 = c^2;
cada1f84 = 2*cada1f83;
cada1f85dx = cada1f82dx./cada1f84;
cada1f85 = cada1f82/cada1f84;
cada1tf1 = cada1f85(Gator1Data.Index24);
cada1f86dx = exp(cada1tf1(:)).*cada1f85dx;
cada1f86 = exp(cada1f85);
cada1f87dx = cada1f74.*cada1f86dx;
cada1f87 = cada1f74*cada1f86;
cada1td1 = zeros(8192,1);
cada1td1(Gator1Data.Index25) = cada1f66dx;
cada1td1(Gator1Data.Index26) = cada1td1(Gator1Data.Index26) + cada1f87dx;
cada1f88dx = cada1td1;
cada1f88 = cada1f66 + cada1f87;
cada1f89 = w(5);
cada1f90 = cada1f89*1;
cada1f91 = cada1f90/one;
cada1f92 = cada1f91*1;
cada1f93 = abs(c);
cada1f94 = cada1f93^2;
cada1f95 = Gator1Data.Data5*cada1f94;
cada1f96 = cada1f92/cada1f95;
cada1f97dx = x.dx(5);
cada1f97 = x.f(5);
cada1tempdx = cada1f97dx(Gator1Data.Index27);
cada1f98dx = -cada1tempdx;
cada1f98 = Xt - cada1f97;
cada1f99dx = 2.*cada1f98(:).^(2-1).*cada1f98dx;
cada1f99 = cada1f98.^2;
cada1f100dx = x.dx(11);
cada1f100 = x.f(11);
cada1tempdx = cada1f100dx(Gator1Data.Index28);
cada1f101dx = -cada1tempdx;
cada1f101 = Yt - cada1f100;
cada1f102dx = 2.*cada1f101(:).^(2-1).*cada1f101dx;
cada1f102 = cada1f101.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index29) = cada1f99dx;
cada1td1(Gator1Data.Index30) = cada1td1(Gator1Data.Index30) + cada1f102dx;
cada1f103dx = cada1td1;
cada1f103 = cada1f99 + cada1f102;
cada1f104dx = -cada1f103dx;
cada1f104 = uminus(cada1f103);
cada1f105 = c^2;
cada1f106 = 2*cada1f105;
cada1f107dx = cada1f104dx./cada1f106;
cada1f107 = cada1f104/cada1f106;
cada1tf1 = cada1f107(Gator1Data.Index31);
cada1f108dx = exp(cada1tf1(:)).*cada1f107dx;
cada1f108 = exp(cada1f107);
cada1f109dx = cada1f96.*cada1f108dx;
cada1f109 = cada1f96*cada1f108;
cada1td1 = zeros(10240,1);
cada1td1(Gator1Data.Index32) = cada1f88dx;
cada1td1(Gator1Data.Index33) = cada1td1(Gator1Data.Index33) + cada1f109dx;
cada1f110dx = cada1td1;
cada1f110 = cada1f88 + cada1f109;
cada1f111 = w(6);
cada1f112 = cada1f111*1;
cada1f113 = cada1f112/one;
cada1f114 = cada1f113*1;
cada1f115 = abs(c);
cada1f116 = cada1f115^2;
cada1f117 = Gator1Data.Data6*cada1f116;
cada1f118 = cada1f114/cada1f117;
cada1f119dx = x.dx(6);
cada1f119 = x.f(6);
cada1tempdx = cada1f119dx(Gator1Data.Index34);
cada1f120dx = -cada1tempdx;
cada1f120 = Xt - cada1f119;
cada1f121dx = 2.*cada1f120(:).^(2-1).*cada1f120dx;
cada1f121 = cada1f120.^2;
cada1f122dx = x.dx(12);
cada1f122 = x.f(12);
cada1tempdx = cada1f122dx(Gator1Data.Index35);
cada1f123dx = -cada1tempdx;
cada1f123 = Yt - cada1f122;
cada1f124dx = 2.*cada1f123(:).^(2-1).*cada1f123dx;
cada1f124 = cada1f123.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index36) = cada1f121dx;
cada1td1(Gator1Data.Index37) = cada1td1(Gator1Data.Index37) + cada1f124dx;
cada1f125dx = cada1td1;
cada1f125 = cada1f121 + cada1f124;
cada1f126dx = -cada1f125dx;
cada1f126 = uminus(cada1f125);
cada1f127 = c^2;
cada1f128 = 2*cada1f127;
cada1f129dx = cada1f126dx./cada1f128;
cada1f129 = cada1f126/cada1f128;
cada1tf1 = cada1f129(Gator1Data.Index38);
cada1f130dx = exp(cada1tf1(:)).*cada1f129dx;
cada1f130 = exp(cada1f129);
cada1f131dx = cada1f118.*cada1f130dx;
cada1f131 = cada1f118*cada1f130;
cada1td1 = zeros(12288,1);
cada1td1(Gator1Data.Index39) = cada1f110dx;
cada1td1(Gator1Data.Index40) = cada1td1(Gator1Data.Index40) + cada1f131dx;
cada1f132dx = cada1td1;
cada1f132 = cada1f110 + cada1f131;
cada1f133dx = -cada1f132dx;
cada1f133 = uminus(cada1f132);
cada1f134 = w(1);
cada1f135 = cada1f134*1;
cada1f136 = cada1f135/one;
cada1f137 = cada1f136*1;
cada1f138 = abs(c);
cada1f139 = cada1f138^2;
cada1f140 = Gator1Data.Data7*cada1f139;
cada1f141 = cada1f137/cada1f140;
cada1f142dx = x.dx(1);
cada1f142 = x.f(1);
cada1tempdx = cada1f142dx(Gator1Data.Index41);
cada1f143dx = -cada1tempdx;
cada1f143 = Xt - cada1f142;
cada1f144dx = 2.*cada1f143(:).^(2-1).*cada1f143dx;
cada1f144 = cada1f143.^2;
cada1f145dx = x.dx(7);
cada1f145 = x.f(7);
cada1tempdx = cada1f145dx(Gator1Data.Index42);
cada1f146dx = -cada1tempdx;
cada1f146 = Yt - cada1f145;
cada1f147dx = 2.*cada1f146(:).^(2-1).*cada1f146dx;
cada1f147 = cada1f146.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index43) = cada1f144dx;
cada1td1(Gator1Data.Index44) = cada1td1(Gator1Data.Index44) + cada1f147dx;
cada1f148dx = cada1td1;
cada1f148 = cada1f144 + cada1f147;
cada1f149dx = -cada1f148dx;
cada1f149 = uminus(cada1f148);
cada1f150 = c^2;
cada1f151 = 2*cada1f150;
cada1f152dx = cada1f149dx./cada1f151;
cada1f152 = cada1f149/cada1f151;
cada1tf1 = cada1f152(Gator1Data.Index45);
cada1f153dx = exp(cada1tf1(:)).*cada1f152dx;
cada1f153 = exp(cada1f152);
cada1f154dx = cada1f141.*cada1f153dx;
cada1f154 = cada1f141*cada1f153;
cada1f155dx = cada1f154dx;
cada1f155 = baseline + cada1f154;
cada1f156 = w(2);
cada1f157 = cada1f156*1;
cada1f158 = cada1f157/one;
cada1f159 = cada1f158*1;
cada1f160 = abs(c);
cada1f161 = cada1f160^2;
cada1f162 = Gator1Data.Data8*cada1f161;
cada1f163 = cada1f159/cada1f162;
cada1f164dx = x.dx(2);
cada1f164 = x.f(2);
cada1tempdx = cada1f164dx(Gator1Data.Index46);
cada1f165dx = -cada1tempdx;
cada1f165 = Xt - cada1f164;
cada1f166dx = 2.*cada1f165(:).^(2-1).*cada1f165dx;
cada1f166 = cada1f165.^2;
cada1f167dx = x.dx(8);
cada1f167 = x.f(8);
cada1tempdx = cada1f167dx(Gator1Data.Index47);
cada1f168dx = -cada1tempdx;
cada1f168 = Yt - cada1f167;
cada1f169dx = 2.*cada1f168(:).^(2-1).*cada1f168dx;
cada1f169 = cada1f168.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index48) = cada1f166dx;
cada1td1(Gator1Data.Index49) = cada1td1(Gator1Data.Index49) + cada1f169dx;
cada1f170dx = cada1td1;
cada1f170 = cada1f166 + cada1f169;
cada1f171dx = -cada1f170dx;
cada1f171 = uminus(cada1f170);
cada1f172 = c^2;
cada1f173 = 2*cada1f172;
cada1f174dx = cada1f171dx./cada1f173;
cada1f174 = cada1f171/cada1f173;
cada1tf1 = cada1f174(Gator1Data.Index50);
cada1f175dx = exp(cada1tf1(:)).*cada1f174dx;
cada1f175 = exp(cada1f174);
cada1f176dx = cada1f163.*cada1f175dx;
cada1f176 = cada1f163*cada1f175;
cada1td1 = zeros(4096,1);
cada1td1(Gator1Data.Index51) = cada1f155dx;
cada1td1(Gator1Data.Index52) = cada1td1(Gator1Data.Index52) + cada1f176dx;
cada1f177dx = cada1td1;
cada1f177 = cada1f155 + cada1f176;
cada1f178 = w(3);
cada1f179 = cada1f178*1;
cada1f180 = cada1f179/one;
cada1f181 = cada1f180*1;
cada1f182 = abs(c);
cada1f183 = cada1f182^2;
cada1f184 = Gator1Data.Data9*cada1f183;
cada1f185 = cada1f181/cada1f184;
cada1f186dx = x.dx(3);
cada1f186 = x.f(3);
cada1tempdx = cada1f186dx(Gator1Data.Index53);
cada1f187dx = -cada1tempdx;
cada1f187 = Xt - cada1f186;
cada1f188dx = 2.*cada1f187(:).^(2-1).*cada1f187dx;
cada1f188 = cada1f187.^2;
cada1f189dx = x.dx(9);
cada1f189 = x.f(9);
cada1tempdx = cada1f189dx(Gator1Data.Index54);
cada1f190dx = -cada1tempdx;
cada1f190 = Yt - cada1f189;
cada1f191dx = 2.*cada1f190(:).^(2-1).*cada1f190dx;
cada1f191 = cada1f190.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index55) = cada1f188dx;
cada1td1(Gator1Data.Index56) = cada1td1(Gator1Data.Index56) + cada1f191dx;
cada1f192dx = cada1td1;
cada1f192 = cada1f188 + cada1f191;
cada1f193dx = -cada1f192dx;
cada1f193 = uminus(cada1f192);
cada1f194 = c^2;
cada1f195 = 2*cada1f194;
cada1f196dx = cada1f193dx./cada1f195;
cada1f196 = cada1f193/cada1f195;
cada1tf1 = cada1f196(Gator1Data.Index57);
cada1f197dx = exp(cada1tf1(:)).*cada1f196dx;
cada1f197 = exp(cada1f196);
cada1f198dx = cada1f185.*cada1f197dx;
cada1f198 = cada1f185*cada1f197;
cada1td1 = zeros(6144,1);
cada1td1(Gator1Data.Index58) = cada1f177dx;
cada1td1(Gator1Data.Index59) = cada1td1(Gator1Data.Index59) + cada1f198dx;
cada1f199dx = cada1td1;
cada1f199 = cada1f177 + cada1f198;
cada1f200 = w(4);
cada1f201 = cada1f200*1;
cada1f202 = cada1f201/one;
cada1f203 = cada1f202*1;
cada1f204 = abs(c);
cada1f205 = cada1f204^2;
cada1f206 = Gator1Data.Data10*cada1f205;
cada1f207 = cada1f203/cada1f206;
cada1f208dx = x.dx(4);
cada1f208 = x.f(4);
cada1tempdx = cada1f208dx(Gator1Data.Index60);
cada1f209dx = -cada1tempdx;
cada1f209 = Xt - cada1f208;
cada1f210dx = 2.*cada1f209(:).^(2-1).*cada1f209dx;
cada1f210 = cada1f209.^2;
cada1f211dx = x.dx(10);
cada1f211 = x.f(10);
cada1tempdx = cada1f211dx(Gator1Data.Index61);
cada1f212dx = -cada1tempdx;
cada1f212 = Yt - cada1f211;
cada1f213dx = 2.*cada1f212(:).^(2-1).*cada1f212dx;
cada1f213 = cada1f212.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index62) = cada1f210dx;
cada1td1(Gator1Data.Index63) = cada1td1(Gator1Data.Index63) + cada1f213dx;
cada1f214dx = cada1td1;
cada1f214 = cada1f210 + cada1f213;
cada1f215dx = -cada1f214dx;
cada1f215 = uminus(cada1f214);
cada1f216 = c^2;
cada1f217 = 2*cada1f216;
cada1f218dx = cada1f215dx./cada1f217;
cada1f218 = cada1f215/cada1f217;
cada1tf1 = cada1f218(Gator1Data.Index64);
cada1f219dx = exp(cada1tf1(:)).*cada1f218dx;
cada1f219 = exp(cada1f218);
cada1f220dx = cada1f207.*cada1f219dx;
cada1f220 = cada1f207*cada1f219;
cada1td1 = zeros(8192,1);
cada1td1(Gator1Data.Index65) = cada1f199dx;
cada1td1(Gator1Data.Index66) = cada1td1(Gator1Data.Index66) + cada1f220dx;
cada1f221dx = cada1td1;
cada1f221 = cada1f199 + cada1f220;
cada1f222 = w(5);
cada1f223 = cada1f222*1;
cada1f224 = cada1f223/one;
cada1f225 = cada1f224*1;
cada1f226 = abs(c);
cada1f227 = cada1f226^2;
cada1f228 = Gator1Data.Data11*cada1f227;
cada1f229 = cada1f225/cada1f228;
cada1f230dx = x.dx(5);
cada1f230 = x.f(5);
cada1tempdx = cada1f230dx(Gator1Data.Index67);
cada1f231dx = -cada1tempdx;
cada1f231 = Xt - cada1f230;
cada1f232dx = 2.*cada1f231(:).^(2-1).*cada1f231dx;
cada1f232 = cada1f231.^2;
cada1f233dx = x.dx(11);
cada1f233 = x.f(11);
cada1tempdx = cada1f233dx(Gator1Data.Index68);
cada1f234dx = -cada1tempdx;
cada1f234 = Yt - cada1f233;
cada1f235dx = 2.*cada1f234(:).^(2-1).*cada1f234dx;
cada1f235 = cada1f234.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index69) = cada1f232dx;
cada1td1(Gator1Data.Index70) = cada1td1(Gator1Data.Index70) + cada1f235dx;
cada1f236dx = cada1td1;
cada1f236 = cada1f232 + cada1f235;
cada1f237dx = -cada1f236dx;
cada1f237 = uminus(cada1f236);
cada1f238 = c^2;
cada1f239 = 2*cada1f238;
cada1f240dx = cada1f237dx./cada1f239;
cada1f240 = cada1f237/cada1f239;
cada1tf1 = cada1f240(Gator1Data.Index71);
cada1f241dx = exp(cada1tf1(:)).*cada1f240dx;
cada1f241 = exp(cada1f240);
cada1f242dx = cada1f229.*cada1f241dx;
cada1f242 = cada1f229*cada1f241;
cada1td1 = zeros(10240,1);
cada1td1(Gator1Data.Index72) = cada1f221dx;
cada1td1(Gator1Data.Index73) = cada1td1(Gator1Data.Index73) + cada1f242dx;
cada1f243dx = cada1td1;
cada1f243 = cada1f221 + cada1f242;
cada1f244 = w(6);
cada1f245 = cada1f244*1;
cada1f246 = cada1f245/one;
cada1f247 = cada1f246*1;
cada1f248 = abs(c);
cada1f249 = cada1f248^2;
cada1f250 = Gator1Data.Data12*cada1f249;
cada1f251 = cada1f247/cada1f250;
cada1f252dx = x.dx(6);
cada1f252 = x.f(6);
cada1tempdx = cada1f252dx(Gator1Data.Index74);
cada1f253dx = -cada1tempdx;
cada1f253 = Xt - cada1f252;
cada1f254dx = 2.*cada1f253(:).^(2-1).*cada1f253dx;
cada1f254 = cada1f253.^2;
cada1f255dx = x.dx(12);
cada1f255 = x.f(12);
cada1tempdx = cada1f255dx(Gator1Data.Index75);
cada1f256dx = -cada1tempdx;
cada1f256 = Yt - cada1f255;
cada1f257dx = 2.*cada1f256(:).^(2-1).*cada1f256dx;
cada1f257 = cada1f256.^2;
cada1td1 = zeros(2048,1);
cada1td1(Gator1Data.Index76) = cada1f254dx;
cada1td1(Gator1Data.Index77) = cada1td1(Gator1Data.Index77) + cada1f257dx;
cada1f258dx = cada1td1;
cada1f258 = cada1f254 + cada1f257;
cada1f259dx = -cada1f258dx;
cada1f259 = uminus(cada1f258);
cada1f260 = c^2;
cada1f261 = 2*cada1f260;
cada1f262dx = cada1f259dx./cada1f261;
cada1f262 = cada1f259/cada1f261;
cada1tf1 = cada1f262(Gator1Data.Index78);
cada1f263dx = exp(cada1tf1(:)).*cada1f262dx;
cada1f263 = exp(cada1f262);
cada1f264dx = cada1f251.*cada1f263dx;
cada1f264 = cada1f251*cada1f263;
cada1td1 = zeros(12288,1);
cada1td1(Gator1Data.Index79) = cada1f243dx;
cada1td1(Gator1Data.Index80) = cada1td1(Gator1Data.Index80) + cada1f264dx;
cada1f265dx = cada1td1;
cada1f265 = cada1f243 + cada1f264;
cada1tf1 = cada1f265(Gator1Data.Index81);
cada1f266dx = 1./cada1tf1(:).*cada1f265dx;
cada1f266 = log(cada1f265);
cada1tf1 = y(Gator1Data.Index82);
cada1f267dx = cada1tf1(:).*cada1f266dx;
cada1f267 = y.*cada1f266;
cada1td1 = cada1f133dx;
cada1td1 = cada1td1 + cada1f267dx;
cada1f268dx = cada1td1;
cada1f268 = cada1f133 + cada1f267;
cada1td1 = zeros(1024,12);
cada1td1(Gator1Data.Index83) = cada1f268dx;
cada1td1 = sum(cada1td1,1);
cada1f269dx = cada1td1(:);
cada1f269 = sum(cada1f268);
out.dx = -cada1f269dx;
out.f = uminus(cada1f269);
%User Line: out = -sum(-(baseline +w(1)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(1)).^2 + ([Yt]-x(7)).^2)/(2*(c.^2)))+w(2)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(2)).^2 + ([Yt]-x(8)).^2)/(2*(c.^2)))+w(3)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(3)).^2 + ([Yt]-x(9)).^2)/(2*(c.^2)))+w(4)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(4)).^2 + ([Yt]-x(10)).^2)/(2*(c.^2)))+w(5)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(5)).^2 + ([Yt]-x(11)).^2)/(2*(c.^2)))+w(6)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(6)).^2 + ([Yt]-x(12)).^2)/(2*(c.^2))))+[y].* log(baseline +w(1)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(1)).^2 + ([Yt]-x(7)).^2)/(2*(c.^2)))+w(2)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(2)).^2 + ([Yt]-x(8)).^2)/(2*(c.^2)))+w(3)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(3)).^2 + ([Yt]-x(9)).^2)/(2*(c.^2)))+w(4)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(4)).^2 + ([Yt]-x(10)).^2)/(2*(c.^2)))+w(5)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(5)).^2 + ([Yt]-x(11)).^2)/(2*(c.^2)))+w(6)*1/one *  1/(2*pi*(abs(c).^2))*exp(-(([Xt]-x(6)).^2 + ([Yt]-x(12)).^2)/(2*(c.^2)))));
out.dx_size = 12;
out.dx_location = Gator1Data.Index84;
end


function ADiGator_LoadData()
global ADiGator_negLikelihood_6_ADiGatorGrd
ADiGator_negLikelihood_6_ADiGatorGrd = load('negLikelihood_6_ADiGatorGrd.mat');
return
end