#  =============== 1. The best 100 colors ==================
#' Select Colors from Predefined Palettes
#'
#' This function selects colors from predefined palettes and allows linear interpolation of colors.
#'
#' Cite from: https://mp.weixin.qq.com/s/4ztyLaLILZ5CJMp94MkeBA
#' @param palette_num An integer from 1 to 100 indicating the chosen palette.
#' @param n An integer specifying the number of colors to generate. If NULL, the original palette is returned.
#' @return A character vector of colors in hexadecimal format.
#' @examples
#' Best100(2)
#' Best100(50, 256)
#' @export
Best100 <- function(palette_num, n = NULL) {
  if (palette_num < 1 || palette_num > 100) {
    stop("Invalid palette number. Please choose a number between 1 and 100.")
  }

  palettes <- list(
    "Scheme1" = c('#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0','#F0027F','#BF5B17','#666666'),
    "Scheme2" = c('#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#666666'),
    "Scheme3" = c('#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928'),
    "Scheme4" = c('#FBB4AE','#B3CDE3','#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#F2F2F2'),
    "Scheme5" = c('#B3E2CD','#FDCDAC','#CBD5E8','#F4CAE4','#E6F5C9','#FFF2AE','#F1E2CC','#CCCCCC'),
    "Scheme6" = c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999'),
    "Scheme7" = c('#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3'),
    "Scheme8" = c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F'),
    "Scheme9" = c('#3182BD','#6BAED6','#9ECAE1','#C6DBEF','#E6550D','#FD8D3C','#FDAE6B','#FDD0A2','#31A354','#74C476','#A1D99B','#C7E9C0','#756BB1','#9E9AC8','#BCBDDC','#DADAEB','#636363','#969696','#BDBDBD','#D9D9D9'),
    "Scheme10" = c('#393B79','#5254A3','#6B6ECF','#9C9EDE','#637939','#8CA252','#B5CF6B','#CEDB9C','#8C6D31','#BD9E39','#E7BA52','#E7CB94','#843C39','#AD494A','#D6616B','#E7969C','#7B4173','#A55194','#CE6DBD','#DE9ED6'),
    "Scheme11" = c('#1F77B4','#AEC7E8','#FF7F0E','#FFBB78','#2CA02C','#98DF8A','#D62728','#FF9896','#9467BD','#C5B0D5','#8C564B','#C49C94','#E377C2','#F7B6D2','#7F7F7F','#C7C7C7','#BCBD22','#DBDB8D','#17BECF','#9EDAE5'),
    "Scheme12" = c('#48648E','#3183BE','#9ECAE1','#49642E','#029194','#2ABDBB','#6C479C','#C47DB5','#9F9AC7','#D9D9D9'),
    "Scheme13" = c('#7ABFE2','#A8D5E9','#F2C76C','#FFD19C','#F1D6A6','#E2C37F','#FDAE6B','#A96769','#CA7073','#ED7E8B','#E64759','#B1ADD6','#C9CBE6','#E2E1F0','#5B5C89','#7576A1'),
    "Scheme14" = c('#8DC691','#7E6875','#9E6762','#637A50','#CC7816','#205479','#F48F83','#512C89','#DFCDE3','#C693BE','#3E9992','#876FAB','#ED5A9E','#F397C0','#1D4A9D','#8DB5CE','#64A83D','#C62128','#F6CA15','#F0861D','#B51D8D','#EE4E21'),
    "Scheme15" = c('#3D3D3D','#00743F','#AECE34','#F5DE24','#F9EA72','#F6C3DB','#CE6BA9','#F28F24','#CB232A','#244C8C','#5DB9DD','#BB6D26','#581F0C','#C2C4C6'),
    "Scheme16" = c('#428AC9','#129392','#FFCC4F','#F37E78','#883A96','#B4B7B7'),
    "Scheme17" = c('#2D2D72','#7FD0F1','#C6D4EA','#744C9D','#7CC892','#E31134','#F9BEB8','#F69335','#000000'),
    "Scheme18" = c('#6C74AF','#7BCFF1','#98B1D9','#A186BD','#98D4AB','#EA627A','#FBD2CE','#F9B474'),
    "Scheme19" = c('#990F0F','#B22C2C','#CC5151','#E57E7E','#FFB2B2','#99540F','#B26F2C','#CC8E51','#E5B17E','#FFD8B2','#6B990F','#85B22C','#A3CC51','#C3E57E','#E5FFB2','#0F6B99','#2C85B2','#51A3CC','#7EC3E5','#B2E5FF','#260F99','#422CB2','#6551CC','#8F7EE5','#BFB2FF'),
    "Scheme20" = c('#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148','#B09C85'),
    "Scheme21" = c('#1F77B4','#FF7F0E','#2CA02C','#D62728','#9467BD','#8C564B','#E377C2','#7F7F7F','#BCBD22','#17BECF'),
    "Scheme22" = c('#5050FF','#CE3D32','#749B58','#F0E685','#466983','#BA6338','#5DB1DD','#802268','#6BD76B','#D595A7','#924822','#837B8D','#C75127','#D58F5C','#7A65A5','#E4AF69','#3B1B53','#CDDEB7','#612A79','#AE1F63','#E7C76F','#5A655E','#CC9900','#99CC00','#A9A9A9','#CC9900','#99CC00','#33CC00','#00CC33','#00CC99','#0099CC','#0A47FF','#4775FF','#FFC20A','#FFD147','#990033','#991A00','#996600','#809900','#339900','#00991A','#009966','#008099','#003399','#1A0099','#660099','#990080','#D60047','#FF1463','#00D68F','#14FFB1'),
    "Scheme23" = c('#FED439','#709AE1','#8A9197','#D2AF81','#FD7446','#D5E4A2','#197EC0','#F05C3B','#46732E','#71D0F5','#370335','#075149','#C80813','#91331F','#1A9993','#FD8CC1'),
    "Scheme24" = c('#CC0C00','#5C88DA','#84BD00','#FFCD00','#7C878E','#00B5E2','#00AF66'),
    "Scheme25" = c('#99CCFF','#CCCCFF','#FFD700','#00CCCC','#AFEEEE','#FF99FF','#FF9999','#FFD700','#99FFCC'),
    "Scheme26" = c('#004586','#FF420E','#FFD320','#579D1C','#7E0021','#83CAFF','#314004','#AECF00','#4B1F6F','#FF950E','#C5000B','#0084D1'),
    "Scheme27" = c('#9C9EDE','#7375B5','#4A5584','#CEDB9C','#B5CF6B','#8CA252','#637939','#E7CB94','#E7BA52','#BD9E39','#8C6D31','#E7969C','#D6616B','#AD494A','#843C39','#DE9ED6','#CE6DBD','#A55194','#7B4173'),
    "Scheme28" = c('#4E79A7','#F28E2B','#E15759','#76B7B2','#59A14F','#EDC948','#B07AA1','#FF9DA7','#9C755F','#BAB0AC'),
    "Scheme29" = c('#4E79A7','#A0CBE8','#F28E2B','#FFBE7D','#59A14F','#8CD17D','#B6992D','#F1CE63','#499894','#86BCB6','#E15759','#FF9D9A','#79706E','#BAB0AC','#D37295','#FABFD2','#B07AA1','#D4A6C8','#9D7660','#D7B5A6'),
    "Scheme30" = c('#795F35','#FBC1AF','#EF2C30','#706D6E','#F4D49D','#FFEE00','#FBC41E','#B45241','#5A69B2','#1AB3EF','#DD9785','#938ABB','#FFF4A4','#86817E','#D8D7D7','#F58971','#E56D4C','#F2833F','#28B45D','#FDE09F','#EE1A95','#F6C4E5','#618F75','#9D4076'),
    "Scheme31" = c('#C70E7B','#FC6882','#007BC3','#54BCD1','#EF7C12','#F4B95A','#009F3F','#8FDA04','#AF6125','#F4E3C7','#B25D91','#EFC7E6','#EF7C12','#F4B95A'),
    "Scheme32" = c('#7DCCD3','#4E7147','#BE9C9D','#F7ECD8','#376597','#9888A5','#DBA662'),
    "Scheme33" = c('#990F26','#B33E52','#CC7A88','#E6B8BF','#99600F','#B3823E','#CCAA7A','#E6D2B8','#54990F','#78B33E','#A3CC7A','#CFE6B8','#0F8299','#3E9FB3','#7ABECC','#B8DEE6','#3D0F99','#653EB3','#967ACC','#C7B8E6','#333333','#666666','#999999','#CCCCCC'),
    "Scheme34" = c('#FCE94F','#EDD400','#C4A000','#8AE234','#73D216','#4E9A06','#FCAF3E','#F57900','#CE5C00','#729FCF','#3465A4','#204A87','#AD7FA8','#75507B','#5C3566','#E9B96E','#C17D11','#8F5902','#EF2929','#CC0000','#A40000','#FFFFFF','#EEEEEC','#D3D7CF','#BABDB6','#888A85','#555753','#2E3436','#000000'),
    "Scheme35" = c('#00B8AA','#374649','#FD625E','#F2C811','#5F6B6D','#8AD4EB','#FE9666','#A66999'),
    "Scheme36" = c('#996600','#CC9900','#FFCC00','#FFFF00','#FFFF99','#FFDB9D','#FFCC66','#FF9933','#FF794B','#FF3300','#990000','#333366','#003399','#0066CC','#0083D7','#0099FF','#3E9ADE','#99CCFF','#B4E2FF','#DEFFFF','#FFCCFF','#CCCCFF','#9999FF','#6666CC','#9999CC','#666699','#006600','#009900','#66CC33','#99FF66','#CCF4CC'),
    "Scheme37" = c('#543005','#8C510A','#BF812D','#DFC27D','#F6E8C3','#F5F5F5','#C7EAE5','#80CDC1','#35978F','#01665E','#003C30'),
    "Scheme38" = c('#8E0152','#C51B7D','#DE77AE','#F1B6DA','#FDE0EF','#F7F7F7','#E6F5D0','#B8E186','#7FBC41','#4D9221','#276419'),
    "Scheme39" = c('#40004B','#762A83','#9970AB','#C2A5CF','#E7D4E8','#F7F7F7','#D9F0D3','#A6DBA0','#5AAE61','#1B7837','#00441B'),
    "Scheme40" = c('#7F3B08','#B35806','#E08214','#FDB863','#FEE0B6','#F7F7F7','#D8DAEB','#B2ABD2','#8073AC','#542788','#2D004B'),
    "Scheme41" = c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'),
    "Scheme42" = c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#FFFFFF','#E0E0E0','#BABABA','#878787','#4D4D4D','#1A1A1A'),
    "Scheme43" = c('#A50026','#D73027','#F46D43','#FDAE61','#FEE090','#FFFFBF','#E0F3F8','#ABD9E9','#74ADD1','#4575B4','#313695'),
    "Scheme44" = c('#A50026','#D73027','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#D9EF8B','#A6D96A','#66BD63','#1A9850','#006837'),
    "Scheme45" = c('#9E0142','#D53E4F','#F46D43','#FDAE61','#FEE08B','#FFFFBF','#E6F598','#ABDDA4','#66C2A5','#3288BD','#5E4FA2'),
    "Scheme46" = c('#F7FBFF','#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6','#2171B5','#08519C','#08306B'),
    "Scheme47" = c('#F7FCFD','#E5F5F9','#CCECE6','#99D8C9','#66C2A4','#41AE76','#238B45','#006D2C','#00441B'),
    "Scheme48" = c('#F7FCFD','#E0ECF4','#BFD3E6','#9EBCDA','#8C96C6','#8C6BB1','#88419D','#810F7C','#4D004B'),
    "Scheme49" = c('#F7FCF0','#E0F3DB','#CCEBC5','#A8DDB5','#7BCCC4','#4EB3D3','#2B8CBE','#0868AC','#084081'),
    "Scheme50" = c('#F7FCF5','#E5F5E0','#C7E9C0','#A1D99B','#74C476','#41AB5D','#238B45','#006D2C','#00441B'),
    "Scheme51" = c('#FFFFFF','#F0F0F0','#D9D9D9','#BDBDBD','#969696','#737373','#525252','#252525','#000000'),
    "Scheme52" = c('#FFF5EB','#FEE6CE','#FDD0A2','#FDAE6B','#FD8D3C','#F16913','#D94801','#A63603','#7F2704'),
    "Scheme53" = c('#FFF7EC','#FEE8C8','#FDD49E','#FDBB84','#FC8D59','#EF6548','#D7301F','#B30000','#7F0000'),
    "Scheme54" = c('#FFF7FB','#ECE7F2','#D0D1E6','#A6BDDB','#74A9CF','#3690C0','#0570B0','#045A8D','#023858'),
    "Scheme55" = c('#FFF7FB','#ECE2F0','#D0D1E6','#A6BDDB','#67A9CF','#3690C0','#02818A','#016C59','#014636'),
    "Scheme56" = c('#F7F4F9','#E7E1EF','#D4B9DA','#C994C7','#DF65B0','#E7298A','#CE1256','#980043','#67001F'),
    "Scheme57" = c('#FCFBFD','#EFEDF5','#DADAEB','#BCBDDC','#9E9AC8','#807DBA','#6A51A3','#54278F','#3F007D'),
    "Scheme58" = c('#FFF7F3','#FDE0DD','#FCC5C0','#FA9FB5','#F768A1','#DD3497','#AE017E','#7A0177','#49006A'),
    "Scheme59" = c('#FFF5F0','#FEE0D2','#FCBBA1','#FC9272','#FB6A4A','#EF3B2C','#CB181D','#A50F15','#67000D'),
    "Scheme60" = c('#FFFFE5','#F7FCB9','#D9F0A3','#ADDD8E','#78C679','#41AB5D','#238443','#006837','#004529'),
    "Scheme61" = c('#FFFFD9','#EDF8B1','#C7E9B4','#7FCDBB','#41B6C4','#1D91C0','#225EA8','#253494','#081D58'),
    "Scheme62" = c('#FFFFE5','#FFF7BC','#FEE391','#FEC44F','#FE9929','#EC7014','#CC4C02','#993404','#662506'),
    "Scheme63" = c('#FFFFCC','#FFEDA0','#FED976','#FEB24C','#FD8D3C','#FC4E2A','#E31A1C','#BD0026','#800026'),
    "Scheme64" = c('#0000FF','#3131FF','#6262FF','#9393FF','#C3C3FF','#F6F6FF','#FFD7D7','#FFA5A5','#FF7474','#FF4343','#FF1111'),
    "Scheme65" = c('#00004C','#000091','#0000D6','#2929FF','#8C8CFF','#EFEFFF','#FFACAC','#FF4B4B','#F30000','#C20000','#910000'),
    "Scheme66" = c('#440154','#482475','#404387','#345F8D','#29788E','#20908C','#23A883','#44BF70','#7AD150','#BEDF25','#FDE724'),
    "Scheme67" = c('#0C0786','#41039D','#6B00A8','#900EA3','#B02A8F','#CB4777','#E16560','#F2844B','#FCA635','#FCCF25','#EFF821'),
    "Scheme68" = c('#000003','#170B3A','#410967','#6B176E','#932567','#BB3755','#DD5138','#F37719','#FBA50A','#F6D745','#FCFEA4'),
    "Scheme69" = c('#000003','#150E36','#3A0F6F','#651A80','#8C2980','#B63679','#DE4A67','#F7705B','#FDA06D','#FDCF92','#FBFCBF'),
    "Scheme70" = c('#00224D','#093370','#35456C','#4F576C','#666970','#7C7B78','#948E77','#AFA371','#C8B765','#E4CE52','#FDE737'),
    "Scheme71" = c('#333399','#1177DD','#00B2B3','#32D670','#97EA84','#FCFE98','#CDBF7E','#9A7E63','#987B75','#CBBCB9','#FDFCFC'),
    "Scheme72" = c('#007F00','#00591A','#003233','#000D4D','#001966','#003F7F','#006498','#178BB2','#64B2CC','#B1D8E5','#FCFDFE'),
    "Scheme73" = c('#7F00FF','#4C4FFB','#1996F2','#1ACFE3','#4DF2CE','#80FEB3','#B3F194','#E6CD73','#FF934D','#FF4B26','#FF0000'),
    "Scheme74" = c('#E1D8E2','#A6BFCB','#6D90BF','#5E57B0','#531D7C','#301436','#64194C','#9E3B4F','#C0755D','#D0B39F','#E1D8E1'),
    "Scheme75" = c('#000008','#0000EF','#009EFF','#00EDD0','#00C061','#1EA500','#FFF100','#FF7800','#FF0101','#FFA0A0','#FFFFFF'),
    "Scheme76" = c('#00096B','#0087B7','#03FEFA','#4BE69B','#93CE3A','#CAE51B','#FEFD01','#FFB300','#FF6700','#CD3400','#9D0400'),
    "Scheme77" = c('#000000','#0000CB','#0097FF','#03E2A2','#35C10C','#FAFE01','#FFB600','#FF4000','#BB0037','#7A34A0','#FFFFFF'),
    "Scheme78" = c('#FDFECC','#C9EBB2','#92D8A4','#66C2A4','#52A8A3','#488F9E','#407598','#3E5B93','#41407B','#382D52','#281A2C'),
    "Scheme79" = c('#2539AF','#287FF8','#34BDFF','#6AEAFE','#8BECB0','#CCFDA1','#EFEC79','#FFBE58','#FFA347','#FFBB87','#FFFFFF'),
    "Scheme80" = c('#151D44','#1C4D61','#117D79','#5DA787','#B6CBAF','#FEF6F4','#E7B8A3','#D5796B','#AF4160','#771A5E','#350D36'),
    "Scheme81" = c('#112040','#244397','#2378A3','#4EA8AF','#ACCEC6','#FEFCDA','#D9C661','#90A30C','#34811E','#10552C','#172413'),
    "Scheme82" = c('#040613','#1C1C38','#302F5F','#3D448A','#3E5EA9','#427BB7','#5296C1','#6BB1CB','#8CCBD6','#BBE3E7','#EAFDFD'),
    "Scheme83" = c('#FEEDB0','#FACA8F','#F5A773','#EE845E','#E26253','#CF4456','#B32E5F','#942063','#721A60','#501653','#2F0F3E'),
    "Scheme84" = c('#331418','#531E22','#732724','#90351E','#A54A17','#B66413','#C47F15','#CF9C1D','#D8BA2A','#DEDB3A','#E1FD4B'),
    "Scheme85" = c('#FFFDCD','#EEDF97','#D8C55F','#B8B22E','#8EA20B','#60920C','#32801F','#0F6C2B','#10542C','#193C24','#172313'),
    "Scheme86" = c('#FFF6F4','#DCDFD0','#B6CBAF','#8CB997','#5DA786','#2B937F','#117D79','#18656F','#1C4D61','#1A3652','#151D44'),
    "Scheme87" = c('#042333','#10326B','#40349F','#684496','#8B538D','#B05F82','#D66C6C','#F3824E','#FCA63C','#F7CF45','#E8FA5B'),
    "Scheme88" = c('#E9F6AB','#DBD987','#CFBC66','#C3A14D','#B58740','#A2713B','#8A5E3A','#714E37','#563E30','#3C2F27','#221F1B'),
    "Scheme89" = c('#510979','#9C1BE3','#BF5FF4','#D89EF8','#EACAFC','#FFFFFF','#FFF2BA','#FFDC58','#FFBB33','#FF7C0A','#E6281E'),
    "Scheme90" = c('#001886','#0034F5','#0966FF','#11D6FF','#4AFFB9','#BAFC46','#FED707','#FE6504','#F60101','#850402'),
    "Scheme91" = c('#252B80','#3752A4','#3C6DB4','#48C6EB','#81C998','#BDD638','#FBCD11','#EF5E21','#EB1D22','#7D1315'),
    "Scheme92" = c('#5184B2','#AAD4F8','#F2F5FA','#F1A7B5','#D55276'),
    "Scheme93" = c('#11325D','#365083','#736B9D','#B783AF','#F5A673','#FCDB72'),
    "Scheme94" = c('#020303','#CA4A2E','#E88D2F','#F8F7F2','#3BA595','#20756A'),
    "Scheme95" = c('#46788E','#78B7C9','#FFFFFF','#F6E093','#E58B7B'),
    "Scheme96" = c('#E76254','#EF8A47','#F7AA58','#FFD06F','#FFE6B7','#FFFFFF','#AADCE0','#72BCD5','#528FAD','#376795','#1E466E'),
    "Scheme97" = c('#264654','#277370','#299D91','#8AB07C','#E8C46A','#F3A263','#E56F51'),
    "Scheme98" = c('#01121A','#005F73','#099396','#92D2C1','#EBD7A5','#ED9B00','#CA6702','#BB3F02','#AE2113','#9C2227'),
    "Scheme99" = c('#224767','#6FAEE0','#A1CAED','#D1E6F8','#FEFDFE','#E4CCE4','#CA95C3','#A85EA2','#8E338A'),
    "Scheme100" = c('#27744B','#6C9C7F','#A3BFAD','#D5E1DA','#FBF9FA','#D9C7CC','#B7939D','#935867','#6A1128')

  )

  # 选择对应的配色方案
  chosen_palette <- palettes[[palette_num]]

  if (is.null(n)) {
    # 如果未指定颜色数量，返回原始配色方案
    return(chosen_palette)
  }

  if (n < 1) {
    stop("The number of colors must be at least 1.")
  }

  # 将颜色转换为RGB矩阵
  col2rgb_matrix <- function(col) {
    rgb <- col2rgb(col)
    return(t(rgb) / 255)
  }

  palette_rgb <- col2rgb_matrix(chosen_palette)

  # 生成线性插值的颜色
  generate_interpolated_colors <- function(rgb_matrix, num_colors) {
    interpolate <- function(x, n) {
      approx(seq_along(x), x, n = n)$y
    }

    r <- interpolate(rgb_matrix[, 1], num_colors)
    g <- interpolate(rgb_matrix[, 2], num_colors)
    b <- interpolate(rgb_matrix[, 3], num_colors)

    rgb_colors <- rgb(r, g, b, maxColorValue = 1)
    return(rgb_colors)
  }

  colors <- generate_interpolated_colors(palette_rgb, n)

  return(colors)
}



















#  =============== 2. PlotViolin ==================
#' PlotViolin
#'
#' This function creates a violin plot with optional significance testing.
#'
#' @param df Data frame containing the data to be plotted.
#' @param x Variable name for grouping (column name).
#' @param y Numeric variable name (column name).
#' @param comparisons List of groups to compare. If not provided, no significance testing is performed.
#' @param fill.col Vector of colors for filling (default is colorRampPalette(brewer.pal(9, "Set1"))(6)).
#' @param color Vector of background colors (default is colorRampPalette(brewer.pal(11, "BrBG"))(30)).
#' @param title Title of the plot.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param angle_x_text Angle of the x-axis text (default is 45).
#' @param legend_position Position of the legend (default is "none").
#' @param signif_test Method for significance testing (default is "t.test").
#' @param signif_map Whether to use asterisks to show significance (default is TRUE).
#' @param signif_tip_length Length of the significance markers (default is c(0.01)).
#' @param x_limits Limits for the x-axis (default is NULL, which means ggplot2 will determine the limits), e.g., c(0, 20).
#' @param y_limits Limits for the y-axis (default is NULL, which means ggplot2 will determine the limits).
#' @param x_breaks Interval for the x-axis breaks (default is 1).
#' @param y_breaks Interval for the y-axis breaks (default is 1).
#' @return A ggplot object.
#' @examples
#' df <- data.frame(samples = rep(c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2"), each = 10),
#'                  values = rnorm(60))
#' comparisons <- list(c("A_1", "A_2"), c("B_1", "B_2"), c("C_1", "C_2"))
#' PlotViolin(df, x = "samples", y = "values", comparisons = comparisons)
#' @export
#' @import ggplot2
#' @import ggpubr
#' @import ggsignif
#' @import tidyverse
#' @import ggprism
#' @import vioplot
#' @import RColorBrewer
#' @import grid
#' @import scales
PlotViolin <- function(df, x, y, comparisons = NULL,
                       fill.col = colorRampPalette(brewer.pal(9, "Set1"))(6),
                       color = colorRampPalette(brewer.pal(11, "BrBG"))(30),
                       title = NULL, xlab = NULL, ylab = NULL,
                       angle_x_text = 45, legend_position = "none",
                       signif_test = "t.test", signif_map = TRUE,
                       step_increase=0.1,signif_tip_length = c(0.01),
                       x_limits = NULL, y_limits = NULL,
                       x_breaks = 1, y_breaks = 1) {

  # Create the plot
  p <- ggplot(df, aes_string(x = x, y = y, fill = x)) +
    geom_violin(trim = TRUE, position = position_dodge(width = 0.1), scale = 'width') +
    geom_boxplot(alpha = 1, outlier.size = 0, size = 0.3, width = 0.2, fill = "white") +
    stat_summary(fun = "mean", geom = "point", shape = 21, size = 2, fill = "blue") +
    labs(x = xlab, y = ylab, title = title) +
    theme_prism() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(),
          legend.position = legend_position,
          axis.text = element_text(color = 'black', size = 12),
          legend.text = element_text(color = 'black', size = 12),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = angle_x_text, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = fill.col)

  # Set x and y axis limits if provided
  if (!is.null(x_limits)) {
    p <- p + scale_x_continuous(limits = x_limits, breaks = seq(x_limits[1], x_limits[2], x_breaks))
  }
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits, breaks = seq(y_limits[1], y_limits[2], y_breaks))
  }

  # Add significance layer if comparisons are provided
  if (!is.null(comparisons)) {
    p <- p + geom_signif(comparisons = comparisons,
                         map_signif_level = signif_map,
                         test = signif_test,
                         step_increase = step_increase,
                         tip_length = signif_tip_length,
                         size = 0.8, color = "black")
  }

  return(p)
}

PlotViolin(iris,x = 'Species','Sepal.Length',y_limits = c(5,10),y_breaks = 0.5,comparisons = list(c('setosa','versicolor'),c('versicolor','virginica')),signif_tip_length = 0.0)
