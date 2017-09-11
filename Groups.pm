#!/usr/bin/perl

use Switch;

package Groups;


my @h3_host_shift = (2,	3,	4,	9,	10,	11,	14,	16,	18,	19,	20,	22,	23,	25,	47,	66,	69,	73,	78,	79,	83,	97,	98,	99,	108,	110,	137,	142,	153,	159,	160,	161,	162,	176,	179,	206,	208,	212,	229,	230,	238,	244,	260,	264,	285,	291,	323,	328,	329,	400,	402,	462,	492,	506,	541);
my @h1_host_shift = (2, 9, 14, 15, 22, 47, 61, 62, 71, 73, 78, 85, 88, 89, 97, 100, 102, 113, 130, 132, 138, 144, 146, 149, 151, 153, 154, 155, 157, 167, 168, 169, 171, 172, 173, 176, 177, 180, 186, 200, 201, 202, 203, 205, 206, 209, 210, 211, 212, 218, 223, 227, 230, 232, 235, 238, 240, 243, 251, 252, 257, 261, 268, 274, 275, 277, 278, 285, 286, 288, 289, 293, 294, 299, 302, 314, 323, 324, 325, 326, 331, 337, 389, 415, 420, 434, 453, 459, 466, 470, 515, 541, 542, 548);
my @n1_host_shift = (3, 5, 8, 12, 13, 14, 16, 20, 26, 29, 34, 40, 41, 42, 43, 46, 47, 51, 52, 53, 59, 64, 66, 67, 69, 70, 71, 72, 74, 75, 76, 78, 79, 80, 81, 82, 83, 85, 93, 95, 99, 101, 105, 111, 114, 116, 136, 149, 157, 189, 195, 200, 206, 210, 211, 214, 220, 221, 222, 223, 232, 241, 250, 257, 258, 263, 264, 267, 273, 274, 285, 287, 288, 289, 309, 311, 329, 339, 340, 341, 351, 354, 355, 365, 367, 369, 382, 386, 388, 390, 393, 394, 396, 427, 430, 432, 434, 451, 454, 455);
my @n2_host_shift = (7, 9, 19, 22, 24, 26, 28, 31, 33, 38, 39, 40, 41, 42, 44, 45, 48, 50, 51, 52, 57, 58, 59, 60, 62, 66, 69, 70, 72, 73, 77, 79, 81, 83, 85, 86, 93, 95, 100, 113, 116, 125, 126, 143, 147, 149, 150, 155, 187, 192, 199, 206, 210, 212, 216, 220, 221, 234, 238, 257, 267, 275, 283, 284, 286, 290, 296, 305, 308, 310, 311, 312, 313, 315, 328, 331, 332, 336, 338, 342, 347, 356, 360, 367, 368, 369, 370, 378, 380, 381, 384, 385, 386, 390, 393, 396, 399, 400, 401, 403, 415, 431, 435, 437, 445, 466);
my @h1_host_shift_001 = (203, 168, 299, 251, 288, 201, 167, 252, 302, 62, 9, 238, 314, 324, 275, 285, 154, 172, 176, 459, 420, 2, 211, 202, 130, 470, 274, 257, 14, 323, 89, 294, 261, 235, 100, 286, 415, 200, 206, 15, 85, 78, 210, 71, 453, 466, 337, 22);
my @h3_host_shift_001 = (16, 108, 229, 244, 79, 73, 83, 161, 260, 20, 9 );
my @n2_host_shift_001 = (386, 384, 381, 328, 83, 70, 81, 192, 51, 147, 125, 283, 41, 286, 77, 72, 378, 331, 126, 155, 50, 62, 338, 369, 60, 315, 216, 399, 396) ;
my @n1_host_shift_001 = (189, 382, 214, 340, 311, 274, 157, 430, 455, 74, 341, 220, 221, 288, 351, 80, 264, 289, 365, 339, 52, 46, 59, 42, 34, 47, 393, 427, 3, 67, 309, 329, 29);

# Caton 1982
my @h1_epitopes = qw(141 142 171 173 175 176 178 179 180 169 172 205 206 209 211 153 156 158 182 186 195 220 237 238 253 286 87 88 90 91 92 132);
# Hensley 2009
my @h1_increased_binding = qw(141 142 171 178 179 180 169 193 205 209 211 156 237 87 88 132 257 175 158 106);
my @n1_epitopes = qw(380 381 382 383 384 385 386 388 389 390 393 397 398 199 200 201 202 223 329 330 332 331 333 336 337 339 340 341 343 344 356 363 364 365 366 367);

# Whiley 
my @h3_epitopes = qw( 138 140 142 147 148 146 149 151 153 154 156 158 159 160 161 162 166 168 184 144 145 171 172 173 174 175 176 179 180 181 202 203 204 205 206 208 209 210 212 213 214 60 61 62 63 64 66 67 69 70 289 291 292 294 295 296 310 313 315 316 320 321 323 324 325 326 327 328 112 118 119 133 137 183 186 187 188 189 190 191 192 193 195 198 217 219 223 224 225 228 229 230 231 232 233 234 235 242 243 244 245 246 254 256 258 260 262 263 264 73 75 78 79 83 91 94 96 97 98 99 102 103 104 107 108 110 125 276 277 278 281 );
# Neher Prediction,  dynamics  and visualization of antigenic phenotypes of seasonal influenza viruses, Table S2 (all sites) Group size = 33
my @h3_antigenic_neher = (78, 160, 172, 174, 212, 292, 137, 140, 149, 158, 98, 99, 147, 315, 41, 91, 171, 205, 278, 151, 161, 156, 209, 202, 241, 173, 206, 242, 276, 175, 228, 258, 128);
# Shih 2007
my @h3_shih_epitopes = qw(66 69 70 137 138 140 142 147 149 151 153 158 159 160 161 162 171 172 173 174 175 176 179 180 188 189 190 202 204 205 206 208 209 212 213 217 223 229 233 242 243 258 260 264 291 292 294 315 323 );
my @n2_epitopes = qw(383 384 385 386 387 389 390 391 392 393 394 396 399 400 401 403 197 198 199 200 221 222 328 329 330 331 332 334 336 338 339 341 342 343 344 346 347 357 358 359  366 367 368 369 370);
my @n1_wan_epitopes = qw(248, 249, 250, 273, 309, 338, 339, 341, 343, 396, 397, 456);
#h1 Huang (antigenic), as is in file Tables_main (Huang + 17, from msa) ( 41 H1N1 HA epitope residues (called natural epitope residues) with statistically significant scores)
my @h1_antigenic = qw( 138 144 145 147 150 158 163 142 170 177 200 203 206 207 208 210 211 52 53 60 288 290 291 294 312 327 111 180 222 226 233 239 241 64 71 86 88 90 97 99 284 );


#h1 Ren, Li, Liu (antigenic) - intersection of two methods
my @h1_antigenic_ren = qw(60 71	88	138	142	144	147	158	204	207	210	222	338);

my @h1_antigenic_Huang_and_host_shift = qw(71	200	203	206	210	211	288	294	);

# before oct 2016
#my @n1_surface = (34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 86, 88, 89, 90, 93, 94, 95, 111, 118, 126, 127, 128, 136, 141, 143, 146, 147, 148, 149, 150, 151, 152, 154, 162, 163, 165, 172, 174, 189, 191, 199, 200, 201, 202, 209, 210, 211, 214, 215, 217, 218, 220, 221, 222, 223, 224, 226, 236, 237, 247, 248, 249, 250, 251, 252, 255, 258, 260, 261, 263, 264, 265, 266, 267, 268, 269, 271, 273, 274, 275, 279, 285, 286, 287, 288, 290, 296, 297, 298, 304, 306, 307, 308, 309, 311, 312, 313, 314, 326, 328, 329, 330, 331, 332, 333, 334, 335, , , 336, 337, 338, 340, 341, 349, 350, 351, 352, 360, 361, 362, 363, 364, 365, 372, 374, 375, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 389, 391, 392, 393, 394, 395, 398, 406, 407, 408, 413, 414, 415, 416, 417, 426, 427, 429, 430, 431, 432, 433, 434, 435, 437, 439, 450, 451, 452, 454, 455, 456, 457, 461, 463, 464, 465, 467, 468);
# oct 2016 (espript 2HTY numeration to Krya numeration, surface+semi-surface (blue and dark blue))
my @n1_thick_surface = (34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,88,89,90,93,94,95,111,118,126,127,128,136,141,143,146,147,148,149,150,151,152,154,162,163,165,171,173,188,190,198,199,200,201,208,209,210,213,214,216,217,219,220,221,222,223,225,235,236,246,247,248,249,250,251,254,257,259,260,262,263,264,265,266,267,268,270,272,273,274,278,284,285,286,287,289,295,296,297,303,305,307,308,309,311,312,313,314,326,328,329,330,331,332,334,335,336,337,338,339,340,341,343,344,352,353,354,355,363,364,365,366,367,368,375,377,378,381,382,383,384,385,386,387,388,389,390,392,395,396,397,398,399,402,410,411,412,413,414,415,416,417,426,427,429,430,431,432,433,435,437,438,440,451,452,453,455,456,457,458,462,464,465,466,468,469);
# oct 2016 espript dark blue (residues 1-83 are absent from pdb)
my @n1_surface = (83,84,88,90,93,95,127,147,148,149,150,151,188,200,209,210,217,220,221,235,248,260,262,264,265,266,268,270,272,273,285,287,309,311,313,329,332,335,339,340,341,343,344,354,355,366,378,382,384,385,386,388,389,390,396,397,414,416,431,432,434,450,463,465,467,468);
# oct 2016 (espript 2HTY numeration to Krya numeration)
my @n1_internal = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,85,87,91,92,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,115,116,117,119,120,121,122,123,124,125,129,130,131,132,133,134,135,137,138,139,140,142,144,145,153,155,156,157,158,159,160,161,164,166,167,168,169,170,172,174,175,176,177,178,179,180,181,182,183,184,185,186,187,189,191,192,193,194,195,196,197,202,203,204,205,206,207,211,212,215,218,224,226,227,228,229,230,231,232,233,234,237,238,239,240,241,242,243,244,245,252,253,255,256,258,261,269,271,275,276,277,279,280,281,282,283,288,290,291,292,293,294,298,299,300,301,302,304,306,310,315,316,317,318,319,320,321,322,323,324,325,327,333,342,345,346,347,348,349,350,351,356,357,358,359,360,361,362,369,370,371,372,373,374,376,379,380,391,393,394,400,401,403,404,405,406,407,408,409,418,419,420,421,422,423,424,425,428,434,439,441,442,443,444,445,446,447,448,449,450,454,459,460,461,463,467);
# before oct 2016
#my @n1_internal = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,85,87,91,92,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,112,113,114,115,116,117,119,120,121,122,123,124,125,129,130,131,132,133,134,135,137,138,139,140,142,144,145,153,155,156,157,158,159,160,161,164,166,167,168,169,170,172,174,175,176,177,178,179,180,181,182,183,184,185,186,187,189,191,192,193,194,195,196,197,202,203,204,205,206,207,211,212,215,218,224,226,227,228,229,230,231,232,233,234,237,238,239,240,241,242,243,244,245,252,253,255,256,258,261,269,271,275,276,277,279,280,281,282,283,288,290,291,292,293,294,298,299,300,301,302,304,306,310,315,316,317,318,319,320,321,322,323,324,325,327,333,342,345,346,347,348,349,350,351,356,357,358,359,360,361,362,369,370,371,372,373,374,376,379,380,391,393,394,400,401,403,404,405,406,407,408,409,418,419,420,421,422,423,424,425,428,434,439,441,442,443,444,445,446,447,448,449,450,454,459,460,461,463,467,470,471,472);

# oct 2016 (espript 4gzp no ligand tetramer, dark blue for surface and white for internal)
my @n2_surface = (82,83,88,90,93,99,100,107,108,110,111,113,126,127,137,141,143,144,146,147,150,151,154,155,162,164,169,170,172,173,176, 187,196,197, 199,200,202,208,209,210,211,214,216,220,221,234,247,249,251,258,259,261,263,264,265,267,269,270,283,284,285,308,309,311,313,328,329,331,332,334,336,338,339,342,344,346,347,357,358,369,381,387,390, 392,394,400,401,413,414,415,416,431,432,434,435,437,450,451,452,454,455,463,465,466,469);
my @n2_internal = (85,87,94,96,97,103,105,106,109,114,115,116,117,119,120,121,122,123,124,125,131,132,133,134,135,136,138,139,140,145,148,156,158,159,167,171,177,178,179,180,181,182,183,184,185,186,188,190,191,192,193,194,201,203,205,207,217,223,225,226,227,228,229,230,231,232,233,235,237,238,239,240,241,242,243,244,252,254,255,256,257,260,272,274,275,276,278,279,280,281,282,287,289,290,291,292,293,297,298,299,300,301,302,303,305,307,312,314,316,317,318,319,320,321,322,323,324,325,327,333,340,345,348,349,350,351,352,353,354,355,360,361,362,363,364,365,366,373,374,375,376,377,379,382,395,397,398,404,405,406,407,408,409,410,411,418,420,421,422,423,424,425,426,427,428,429,436,438,439,440,441,442,443,444,445,446,448,460,462,464,467);
# before oct 2016
#my @n2_surface = (34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 88, 89, 90, 91, 92, 93, 95, 107, 110, 111, 112, 113, 118, 125, 126, 127, 128, 130, 141, 143, 146, 147, 149, 150, 151, 152, 153, 154, 160, 161, 162, 169, 171, 173, 187, 189, 196, 197, 198, 199, 200, 208, 209, 210, 212, 215, 216, 218, 219, 220, 221, 222, 224, 234, 236, 244, 245, 246, 247, 248, 249, 250, 251, 253, 258, 259, 261, 262, 263, 264, 265, 267, 268, 269, 270, 271, 273, 277, 283, 284, 285, 286, 292, 295, 296, 304, 306, 307, 308, 309, 310, 311, 312, 313, 315, 326, 328, 329, 330, 331, 332, 334, 336, 337, 338, 339, 341, 342, 343, 344, 346, 347, 356, 357, 358, 359, 366, 367, 368, 369, 370, 371, 378, 380, 381, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 396, 399, 400, 401, 402, 403, 413, 414, 415, 416, 417, 430, 431, 432, 433, 434, 435, 437, 450, 451, 452, 453, 455, 456, 457, 459, 461, 463, 464, 465, 466, 468, 469, 470);
#my @n2_internal = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 87, 94, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 108, 109, 114, 115, 116, 117, 119, 120, 121, 122, 123, 124, 129, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 142, 144, 145, 148, 155, 156, 157, 158, 159, 163, 164, 165, 166, 167, 168, 170, 172, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 188, 190, 191, 192, 193, 194, 195, 201, 202, 203, 204, 205, 206, 207, 211, 213, 214, 217, 223, 225, 226, 227, 228, 229, 230, 231, 232, 233, 235, 237, 238, 239, 240, 241, 242, 243, 252, 254, 255, 256, 257, 260, 266, 272, 274, 275, 276, 278, 279, 280, 281, 282, 287, 288, 289, 290, 291, 293, 294, 297, 298, 299, 300, 301, 302, 303, 305, 314, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 327, 333, 335, 340, 345, 348, 349, 350, 351, 352, 353, 354, 355, 360, 361, 362, 363, 364, 365, 372, 373, 374, 375, 376, 377, 379, 382, 395, 397, 398, 404, 405, 406, 407, 408, 409, 410, 411, 412, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 436, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 454, 458, 460, 462, 467);

#h3 antigentic Steinbruck Table 1
my @h3_antigenic  = qw( 138 160 171 223 161 205 233 294 66 153 174 276 140 151 230 278 78 172 212 292 41 91 99 147 202 218 238 241 19 204 69 180 190 209 217 229 246 98 149 159 162 176 213 18 70 188 260 206 242 );

# h3 antigenic change Koel
my @h3_antigenic_koel = (161, 171, 172, 174, 175, 205, 209);
my @h3_antigenic_smith = (138, 160, 153, 161, 149, 159, 162, 140, 147, 171, 204, 180, 205, 209, 174, 172, 176, 213, 175, 206, 212, 294, 66, 70, 223, 190, 217, 229, 233, 246, 188, 260, 292, 98, 276, 78, 91, 99, 69);
# non-epitopic sites added 
my @h3_antigenic_smith_full = (41,218,238,241, 138, 160, 153, 161, 149, 159, 162, 140, 147, 171, 204, 180, 205, 209, 174, 172, 176, 213, 175, 206, 212, 294, 66, 70, 223, 190, 217, 229, 233, 246, 188, 260, 292, 98, 276, 78, 91, 99, 69);

# h3 "truly antigenic" Steinbruck (subset of Table 1 sites)
my @h3_antigenic_steinbruck  = (138, 160, 171, 223, 161, 205, 233, 294, 66, 153, 174, 276, 140, 151, 278, 230, 78, 172, 212, 292, 41, 91, 99, 147, 202, 218, 238, 241);

#my @h1_pocket_closest = (207,239,166,211,148,111,208,206,234,149,235,196,203,210,151,150,241,204,237,236,240,238,205,209,242,110,197,195,147,202,212,165,152,112,167,233,158,243,244,198,265,159,199,263,201,153,169,262,160);
my @h1_pocket_closest = (207,239,166,111,149,196,203,241,147,168,240,242,110,169,146,197,195,148,206,238,150,204,208,202,167,165,112,243,244,235,198,265,158,199,263,144,205,151,236,159,201,237,152,145,200,143,113,209,245);
my @h3_pocket_closest = (152,241,242,206,238,210,151,114,153,169,165,205,209,42,29,28,211,239,150,207,237,168,113,154,170,115,240,243,212,208,166,43,162,164,161,171,268,175,27,203,204,163,149,155,244,167,172,41,30,148,236,44,156,269);
my @n2_pocket_closest = (151,277,406,150,407,276,278,152,405,226,350,425,291,225,292);
my @n2_bigger_pocket_closest = (292,152,119,276,224,227,406,118,277,371,117,225,407,275,278,372,370,291,223,153,151,228,405,293,226,120,350,425,180,404,181,242,178,133,134,300,241,440,179,441,365,222,150,349);
my @n1_pocket_closest = (402,151,278,152,403,401,277,279,150,227,425,347,292,424,226,228);

my @h1_surface = (28,30,38,39,40,42,52,53,54,56,62,63,64,68,71,83,86,91,92,100,101,103,111,132,136,137,138,141,142,146,148,150,155,156,157,158,171,172,173,176,178,182,184,186,200,201,202,205,206,211,212,221,227,232,235,237,238,252,253,277,278,283,285,287,289,290,292,294,299,303,305,313,324,326,327,339,340,344,350,351,354,358,359,361,362,370,372,373,375,377,381,382,385,386,389,392,396,400,403,404,406,410,412,415,416,425,464,467,470,471,474,476,477,478,484,485,486,488,490,493,497,498,499,501,502,503);
my @h1_internal = (20,22,23,24,25,26,33,35,36,41,43,44,49,50,58,59,66,67,69,72,74,75,76,77,78,79,80,81,82,84,87,93,95,96,97,98,104,105,107,108,109,112,117,118,119,121,122,125,128,133,134,140,143,149,151,152,160,161,163,164,165,166,167,174,177,183,189,190,191,192,193,194,195,196,197,198,199,204,208,213,215,216,217,218,219,226,231,233,241,242,243,245,246,247,248,250,254,256,258,259,260,262,263,264,265,266,267,270,271,272,273,280,282,284,293,295,296,297,298,300,301,302,307,308,309,315,316,319,320,323,325,328,330,331,333,334);

my @h3_surface = (25,37,38,41,43,47,48,49,54,61,62,64,66,71,73,79,94,97,98,99,107,108,110,112,120,120,140,142,144,147,148,149,151,153,156,158,159,160,161,173,174,175,176,178,179,181,183,187,188,189,204,205,206,208,209,214,215,224,228,230,238,240,241,255,256,277,278,279,280,285,287,289,292,294,301,305,307,326,328,329,340,341,352,353,356,361,363,364,372,374,375,376,377,383,384,387,391,394,398,402,403,405,406,414,416,418,427,466,466,480,488,491,492,499,500,501,503,505,506,509,510,513,517,518);
my @h3_internal = (27,29,31,32,33,35,42,44,45,50,52,53,58,59,60,67,68,72,75,77,80,82,83,84,85,86,87,88,89,92,95,100,102,103,104,105,106,113,114,115,118,123,124,125,126,127,128,129,131,132,133,134,136,141,143,146,155,163,164,166,167,168,169,170,177,180,182,186,192,193,194,195,196,197,198,199,200,201,202,207,211,216,218,219,220,221,222,229,231,236,244,245,246,247,248,251,253,257,259,260,261,263,265,266,267,268,269,270,272,273,274,281,282,283,284,286,288,291,297,298,299,302,303,304,310,311,318,319,321,322,325,330,332,333,335,336,338,349,350,351,355,358,359,362,367,368,369,373,385,386,389,393,396,411,420,423,425,426,428,429,430,432,434,435,436,437,438,440,441,444,445,446,448,449,452,453,454,455,456,457,458,459,460,463,464,467,469,471,474,476,477,481,483,485,486,487,489,493,494,497,502,508,511,512,515);


our @h1_leading_kr = (4,11,13,16,52,60,73,74,86,88,90,91,97,99,111,113,128,144,151,156,157,162,169,170,171,172,173,178,182,184,199,202,203,205,206,207,209,210,232,240,261,268,269,283,287,289,290,293,326,361,415,488,489);
our @h1_trailing_kr = (3,6,7,11,52,53,64,89,91,99,101,111,129,137,142,144,148,150,156,157,158,162,165,169,172,178,179,185,186,195,197,199,200,201,203,207,231,236,243,251,253,269,274,278,289,290,293,299,324,331,361,389,390,398,415,422,455,467,470,489,510,514,515,526,562);
#before 23.03.2017
# our @h3_leading_kr = (27,35,57,82,89,94,107,115,138,153,156,163,165,167,169,172,175,176,177,187,188,190,191,192,195,204,218,221,222,224,225,229,234,249,251,254,258,259,293,294,307,308,310,379,393,407,418,482,484,561);
my @h3_leading_kr = (11,19,41,66,73,78,91,99,122,137,140,147,149,151,153,156,159,160,161,171,172,174,175,176,179,188,202,205,206,208,209,213,218,233,235,238,242,243,277,278,291,292,294,363,377,391,402,466,468,545);
#before 23.03.2017
# our @h3_trailing_kr = (19,26,27,29,32,59,65,77,79,80,81,82,83,85,88,89,107,117,120,124,126,128,138,144,160,163,169,170,172,174,182,189,191,192,195,196,203,204,205,206,224,225,226,231,233,234,239,241,246,248,252,254,255,257,258,261,276,280,292,301,305,308,311,323,355,358,393,407,417,418,450,458,482,500,521,532,554,562,579);
my @h3_trailing_kr = (3,10,11,13,16,43,49,61,63,64,65,66,67,69,72,73,91,101,104,108,110,112,122,128,144,147,153,154,156,158,166,173,175,176,179,180,187,188,189,190,208,209,210,215,217,218,223,225,230,232,236,238,239,241,242,245,260,264,276,285,289,292,295,307,339,342,377,391,401,402,434,442,466,484,505,516,538,546,563);
our @n2_leading_kr = (18,20,23,30,52,93,143,150,194,197,199,208,216,220,221,249,265,307,308,310,313,328,336,339,344,346,368,369,370,372,381,385,387,390,432);
our @n2_trailing_kr = (2,4,5,9,27,30,40,44,45,50,56,65,77,82,83,120,127,147,148,149,151,155,210,216,220,238,248,251,258,263,265,269,302,307,309,310,312,328,329,334,335,338,339,342,347,372,386,392,400,402,403,414,416,432,433,434,455,464);
our @n1_leading_kr = (15,17,23,34,45,64,70,78,105,173,200,214,222,234,248,249,250,254,270,274,275,287,329,332,336,339,344,352,354,367,369,382,390,396,418,427,430,434,451);
our @n1_trailing_kr = (15,17,21,23,38,39,40,42,45,47,48,52,57,67,68,70,73,77,81,82,83,93,100,101,114,130,147,149,188,200,249,254,259,262,264,267,270,273,275,329,331,340,346,352,364,366,367,390,416,418,419,427,435,452,453,455,462);


my @h3_leading_neva = (18,176,159,175,243,242,162,260,264,161,377,213,202,19,292,69,208,178,531,138,172);
my @h3_trailing_neva = (97,79,236,342,538,13,179,468,289,264,235,173,400);
my @n2_leading_neva = (126,56,249,332,399,264,431,215,43,267,248,220,290,328,313,46,208,432,172,372,69,308,370);
my @n2_trailing_neva = (401,372,155,335,220,339,432,430,44,127,263);

my @h1_leading_neva = (4,113,171,257,178,138,290,52,415,184,223,16,102,86,199,157,162,361,200,111,74,145,221,169,64);
my @h1_trailing_neva = (91,200,232,169,268);
my @n1_leading_neva = (163,263,388,6,149,59,78,14,80,101,427,386,200);
my @n1_trailing_neva = (434,275,15,267,83,287);


my @h3_deps_evolving = qw(10 61 151 161 171 174 245 264 347);
## 90 240 277 179 - only in egg-adapted before 1979
my @h1_jianpeng_evolving = qw(156 169 171 203 206 210 238 90 240 277 179);

my @h1_wenfu_evolving = qw(98 110 157 178 202 203 238 176);

#


sub get_predefined_groups_and_names_for_protein {
	my $prot = shift;
	my $length = shift;
	my @groups;
	my @names;
	my @real_groups_and_names = only_groups($prot);
	my @sites = (1..$length);
	my @groups_and_names = prepare_groups_and_names($real_groups_and_names[0], $real_groups_and_names[1], \@sites);
	return @groups_and_names;
}


## sites with median stat p-values < 0.01 (from Groups_and_sites.xls), restriction = 50
my @h1_nsyn_epi = (171,184,202,287,4,90,73,151,88,179);
my @h1_nsyn_env = (206,289,13,144);
my @h1_nsyn_epi_005 = (171,184,202,361,287,4,90,73,151,88,179,201,415,238,232,269); # p-values < 0.05 
my @h1_nsyn_env_005 = (155,209,470,169,97,223,207,169,12,206,289,13,144);
my @h1_syn_epi_005 = (484,108,149,478,176,242,453,453,98,248,10,59,170,354,123,561,60,280,525,379,254,543,89,2,539,514,476,40,241,55,61,495);
my @h1_syn_env_005 = (304,495,51,250,451,49,507,498,70,356,290,352,482,156,58,249,506,134,550,501,170,478,228);
my @h1_syn_epi = (354,123,561,60,280,525,379,254,543,89,2,539,514,476,40,241,55,61,495);
my @h1_syn_env = (304,495,51,250,451,49,507);


sub best_sites {
	my $prot = shift;
	my $state = shift;
	my @array;
	switch ($prot) {
		case "h1" {	
			if ($state == "nsyn"){
				push @array, @h1_nsyn_env;
				push @array, @h1_nsyn_epi;
			}
			if ($state == "syn"){
				push @array, @h1_syn_env;
				push @array, @h1_syn_epi;
			}
					  }
		else {
			print "Panic! no best sites defined for protein $prot\n";
			die;
		}			  
	}
	return @array;
}



sub get_fake_predefined_groups_and_names_for_protein {
	my $prot = shift;
	my $length = shift;
#	my $exclude = shift;
	my $state = shift;
	my @groups;
	my @real_groups_and_names = only_groups($prot);
	my @allsites = (1..$length);
	my @sites;
#	if ($exclude){
#		my @sites_to_exclude = best_sites($prot, $state);
#		my %exclude_map = map {$_ => 1} @sites_to_exclude;
#		@sites = grep {not $exclude_map{$_}} @allsites;
#	}

#	else {
		@sites = @allsites;
#	}
	my @heap_of_sites = @sites; # it's a copy, yes
	foreach my $g(@{$real_groups_and_names[0]}){
		my $group_size = scalar @{$g};
		my @group;
		
  		for ( 1..$group_size ){
  			push @group, splice @heap_of_sites, rand @heap_of_sites, 1;
  		}
  		print "\n new group ";
  		foreach my $s (@group){
  			print $s."\t";
  		}
  		
  		push @groups, \@group;
	}
	my @groups_and_names = prepare_groups_and_names(\@groups, $real_groups_and_names[1], \@sites);
	return @groups_and_names;
	
}

# without complement
sub only_groups_debugging {
	my $prot = shift;
	my @groups;
	my @names;
	if ($prot eq "h1"){
		@groups = ( \@h1_internal);
		@names = ( "internal");
	}
	return (\@groups, \@names)
}

sub only_groups {
	my $prot = shift;
	my @groups;
	my @names;
	if ($prot eq "h1"){
		@groups = (\@h1_increased_binding, \@h1_antigenic, \@h1_pocket_closest, \@h1_surface, \@h1_internal, \@h1_host_shift_001, \@h1_leading_kr, \@h1_trailing_kr, \@h1_antigenic_ren);
		@names = ("increased_binding", "antigenic", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr", "antigenic_ren");
	}
	elsif ($prot eq "h3"){
		@groups = (\@h3_shih_epitopes, \@h3_antigenic, \@h3_antigenic_smith, \@h3_antigenic_smith_full, \@h3_antigenic_koel, \@h3_antigenic_neher, \@h3_pocket_closest, \@h3_surface, \@h3_internal, \@h3_host_shift_001, \@h3_leading_kr, \@h3_trailing_kr);
		@names = ("shih_epitopes", "antigenic", "antigenic_smith", "antigenic_smith_full", "antigenic_koel", "antigenic_neher", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr");	
	}
	elsif ($prot eq "n1"){
		@groups = (\@n1_epitopes, \@n1_pocket_closest, \@n1_surface, \@n1_thick_surface, \@n1_internal, \@n1_host_shift_001, \@n1_leading_kr, \@n1_trailing_kr);
		@names = ("epitopes", "pocket_closest", "surface", "thick_surface","internal", "host_shift_001", "leading_kr", "trailing_kr");
	}
	elsif ($prot eq "n2"){
		@groups = (\@n2_epitopes, \@n2_pocket_closest, \@n2_surface, \@n2_internal, \@n2_host_shift_001, \@n2_leading_kr, \@n2_trailing_kr, \@n2_decreasing, \@n2_increasing);
		@names = ("epitopes", "pocket_closest", "surface", "internal", "host_shift_001", "leading_kr", "trailing_kr", "decreasing", "increasing");
	}
	return (\@groups, \@names)
}


sub get_no_groups_for_protein {
	my $prot = shift;
	my $length = shift;
	my @all_sites = (1..$length);
	my @groups;
	my @group_names;
	push @groups, \@all_sites;
	push @group_names, "all";
	return (\@groups, \@group_names);
}

	
sub describe_site{
	my $prot = shift;
	my $site = shift;
	my  @groups_and_names = get_predefined_groups_and_names_for_protein($prot);
	for (my $i = 0; $i < scalar @{$groups_and_names[0]}; $i++){
		if (grep /^$site$/, @{$groups_and_names[0][$i]} ){
			if (! grep /complement$/, $groups_and_names[1][$i] ){
				if ( ! grep /^all$/, $groups_and_names[1][$i]){
				print $groups_and_names[1][$i]." ";
				}
			}
		}
	}
	
}	




sub prepare_groups_and_names {
	my @pregroups = @{$_[0]};
	my @pregroup_names = @{$_[1]};
	my @sites = @{$_[2]};
	my @group_names;
	
	foreach my $group_name(@pregroup_names){
		push @group_names, $group_name;
		push @group_names, $group_name."_complement";
	}
	
	my @groups = prepare_groups(\@pregroups, \@sites);

	#my @all_variable_sites;
	#my %nodes_with_sub = $realdata->{"nodes_with_sub"};
	#foreach my $s(1..565){
	#	if ($nodes_with_sub{$s}){
	#		print " pushed $s\n";
	#		push @all_variable_sites, $s;
	#	}
	#}
	
#	print ("debugging prep_groups length is $length\n");
	push @groups, \@sites;
	push @group_names, "all";
	
	return (\@groups, \@group_names);
}



sub prepare_groups {
	my @groups = @{$_[0]};
	my @sites = @{$_[1]};
	my @final_groups;
	foreach my $pregroup (@groups){
		my %ghash;
		my @group;
		foreach my $s(@{$pregroup}){
			$ghash{$s} = 1;
			push @group, $s;
		}
		my @complement;
		foreach my $s(@sites){
			if (!$ghash{$s}){
				push @complement, $s;
			}
		}
		push @final_groups, \@group;
		push @final_groups, \@complement;
	}
	return @final_groups;
}	

#16.02 complement to leading is trailing, and vice versa
sub prepare_lt_groups_and_names {
	my $prot = $_[0];
	my @sites = @{$_[1]};
	my @group_names;
	push @group_names, "pure_leading";
	push @group_names, "trailing_as_complement";
	push @group_names, "pure_trailing";
	push @group_names, "leading_as_complement";	

	my @groups = prepare_lt_groups($prot);
	
	push @groups, \@sites;
	push @group_names, "all";
	
	return (\@groups, \@group_names);
}

#16.02 complement to leading is trailing, and vice versa
sub prepare_lt_groups {		
		my $prot = $_[0];
		my @leading;
		my @trailing;
		my @final_groups;
		switch ($prot) {
			case "h1" {	@leading = @h1_leading_kr;
					@trailing = @h1_trailing_kr;
					  }
			case "h3" {	@leading = @h3_leading_kr;
					@trailing = @h3_trailing_kr;
					  }
			case "n1" {	@leading = @n1_leading_kr;
					@trailing = @n1_trailing_kr;
					  }
			case "n2" {	@leading = @n2_leading_kr;
					@trailing = @n2_trailing_kr;
					  }	
			else {	print "Panic in prepare_lt_groups! No such protein $prot!\n";
					die;
				 }			
		}
		my %trailing = map{$_ =>1} @trailing;
		my %leading = map{$_ =>1} @leading;
		
		my @pure_leading = grep(!defined $trailing{$_}, @leading);
		my @pure_trailing = grep(!defined $leading{$_}, @trailing);
	
		push @final_groups, \@pure_leading;
		push @final_groups, \@pure_trailing;
		push @final_groups, \@pure_trailing;
		push @final_groups, \@pure_leading;
		
		return @final_groups;
}

1;