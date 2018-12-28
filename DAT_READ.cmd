/* =========================================================================*
 * OMAPL138_DSP.cmd - Linker Command File for Linking OMAPL138 DSP Programs	*
 * 																			*
 * These linker options are for command line linking only. For IDE linking, *
 * you should set your linker options in Project Properties.				*
 * 		-c					Link Using C Conventions						*
 *		-stack 	0x1000		Software Stack Size								*
 *		-heap	0x1000		Heap Area Size									*
 * =========================================================================*/
-stack 0x100000
// -heap 0x1000
//-e fix_start
/* =========================================================================*
 * 						Specify the System Memory Map						*
 * =========================================================================*/
MEMORY
{
	L2:				o = 0x00800000		l = 0x00040000
	share:			o = 0x80000000		l = 0x00020000	  //128k
	L1P:			o = 0x00E00000		l = 0x00008000
	L1D:			o = 0x00F00000		l = 0x00008000
	DDR2: 	    	o = 0xc2000000  	l = 0x00ffffff    //16M
	DDR2_HYDICE: 	o = 0xc3000000		l = 0x01ffffff
}

/* =========================================================================*
 * 				Specify the Sections Allocation into Memory					*
 * =========================================================================*/
SECTIONS
{
	.cinit		>		DDR2				/* Initialization Tables			*/
	.pinit		>		DDR2				/* C++ Constructor Tables			*/
	.const		>		DDR2				/* Constant Data					*/
	.switch		>		DDR2				/* Jump Tables						*/
	.text:_c_int00 > 0xc2000000
	.text		>		DDR2				/* Executable Code					*/



	.bss		>		DDR2				/* Global & Static Variables 		*/
	.far		>		DDR2				/* Far Global & Static Variables	*/
	.stack		>		DDR2				/* Software System Stack			*/
	.sysmem		>		DDR2				/* Dynamic Memory Allocation Area 	*/

	.cio		>		DDR2				/* C I/O Buffer						*/
	.vecs		>		DDR2				/* Interrupt Vectors				*/


	.MY_MEM:{}  >       DDR2_HYDICE			//导入hex图像的地址
}
