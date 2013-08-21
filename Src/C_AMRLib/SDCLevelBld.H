
#ifndef _SDCLEVELBLD_H_
#define _SDCLEVELBLD_H_

#include <SDCAmr.H>
#include <sdc.h>

extern SDCLevelBld* getSDCLevelBld ();

//
// Builds problem-specific SDCAmrLevels.
//

class SDCLevelBld
{
public:

  virtual ~SDCLevelBld () {};

  virtual sdc_sweeper_t* build (SDCAmr& papa,
                                int     lev,
                                Real    time) = 0;

};

#endif /*_SDCLEVELBLD_H_*/
