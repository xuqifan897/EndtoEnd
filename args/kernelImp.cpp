#include "args.h"

using namespace E2E;
using namespace std;

CCCSkernel::CCCSkernel(int NA): num_angles(NA)
{
    this->angles = new float[NA]{0};
    this->Atheta = new float[NA]{0};
    this->Btheta = new float[NA]{0};
    this->atheta = new float[NA]{0};
    this->btheta = new float[NA]{0};

    this->angles_flag = true;
    this->Atheta_flag = true;
    this->Btheta_flag = true;
    this->atheta_flag = true;
    this->btheta_flag = true;
}

CCCSkernel::~CCCSkernel()
{
    if (this->angles_flag)
        delete this->angles;
    if (this->Atheta_flag)
        delete this->Atheta;
    if (this->Btheta_flag)
        delete this->Btheta;
    if (this->atheta_flag)
        delete this->atheta;
    if (this->btheta_flag)
        delete this->btheta;
}

CCCSkernel::CCCSkernel(CCCSkernel& old): num_angles{old.num_angles}
{
    this->angles_flag = old.angles_flag;
    this->Atheta_flag = old.Atheta_flag;
    this->Btheta_flag = old.Btheta_flag;
    this->atheta_flag = old.atheta_flag;
    this->btheta_flag = old.btheta_flag;

    if (this->angles_flag)
    {
        this->angles = new float[this->num_angles];
        for (int i=0; i<this->num_angles; i++)
            this->angles[i] = old.angles[i];
    }
    if (this->Atheta_flag)
    {
        this->Atheta = new float[this->num_angles];
        for (int i=0; i<this->num_angles; i++)
            this->Atheta[i] = old.Atheta[i];
    }
    if (this->Btheta_flag)
    {
        this->Btheta = new float[this->num_angles];
        for (int i=0; i<this->num_angles; i++)
            this->Btheta[i] = old.Btheta[i];
    }
    if (this->atheta_flag)
    {
        this->atheta = new float[this->num_angles];
        for (int i=0; i<this->num_angles; i++)
            this->atheta[i] = old.atheta[i];
    }
    if (this->btheta_flag)
    {
        this->btheta = new float[this->num_angles];
        for (int i=0; i<this->num_angles; i++)
            this->btheta[i] = old.btheta[i];
    }
}

CCCSkernel::CCCSkernel(CCCSkernel&& old): num_angles(old.num_angles)
{
    this->angles_flag = old.angles_flag;
    this->Atheta_flag = old.Atheta_flag;
    this->Btheta_flag = old.Btheta_flag;
    this->atheta_flag = old.atheta_flag;
    this->btheta_flag = old.btheta_flag;

    if (old.angles_flag)
        this->angles = move(old.angles);
    
    if (old.Atheta_flag)
        this->Atheta = move(old.Atheta);

    if (old.Btheta_flag)
        this->Btheta = move(old.Btheta);

    if (old.atheta_flag)
        this->atheta = move(old.atheta);

    if (old.btheta_flag)
        this->btheta = move(old.btheta);

    old.angles_flag = false;
    old.Atheta_flag = false;
    old.Btheta_flag = false;
    old.atheta_flag = false;
    old.btheta_flag = false;
}

FCBBkernel::FCBBkernel(int ND): num_depths(ND)
{
    this->depths_flag = true;
    this->doses_flag = true;

    this->depths = new float[this->num_depths]{0};
    this->doses = new float[this->num_depths]{0};

    this->A = 0;
    this->B = 0;
    this->a = 0;
    this->b = 0;
}

FCBBkernel::~FCBBkernel()
{
    if (this->depths_flag)
        delete this->depths;
    if (this->doses_flag)
        delete this->doses;
}

FCBBkernel::FCBBkernel(FCBBkernel& old): \
num_depths(old.num_depths), A(old.A), B(old.B), a(old.a), b(old.b)
{
    this->depths_flag = old.depths_flag;
    this->doses_flag = old.doses_flag;

    if (this->depths_flag)
    {
        this->depths = new float[this->num_depths];
        for (int i=0; i<this->num_depths; i++)
            this->depths[i] = old.depths[i];
    }

    if (this->doses_flag)
    {
        this->doses = new float[this->num_depths];
        for (int i=0; i<this->num_depths; i++)
            this->doses[i] = old.doses[i];
    }
}

FCBBkernel::FCBBkernel(FCBBkernel&& old): \
num_depths(old.num_depths), A(old.A), B(old.B), a(old.a), b(old.b)
{
    this->depths_flag = old.depths_flag;
    this->doses_flag = old.doses_flag;

    if (old.depths_flag)
        this->depths = move(old.depths);
    if (old.doses_flag)
        this->doses = move(old.doses);
    
    old.depths_flag = false;
    old.doses_flag = false;
}