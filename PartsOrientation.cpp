/** 
 * file PartsOrientation.cpp
 * author Bruna
 * date Created on 22 de Maio de 2019, 15:29
 * brief Methods to calculate the order to print each part of the object.
 * details During the procedure there are methods to check if there are collisions between the part and the printer in order to create just the sequences that the printer is able to physically print.
 */

#include "PartsOrientation.hpp"

// list of printing sequences
vector<vector<vector<pair<int, int>>>> L_f;


// our printer characteristics
extern printer printer_char;
extern void ANGLES(float vx, float vy, float vz, float *B, float *C, float *lastb, float *lastc);


void my_rotate(float *x, float *y, float angle) { // copy of stl_rotate
  double r;
  double theta;
  double radian_angle;

  radian_angle = (angle / 180.0) * M_PI;

  r = sqrt((*x **x) + (*y **y));
  theta = atan2(*y, *x);
  *x = r * cos(theta + radian_angle);
  *y = r * sin(theta + radian_angle);
}

void
my_stl_rotate(stl_file *stl, float B, float C) {
  int i;
  int j;
  double nx, ny, nz;

  if (stl->error) return;
  
#pragma omp parallel for
  for(i = 0; i < stl->stats.number_of_facets; i++) {
    for(j = 0; j < 3; j++) {
        nx = stl->facet_start[i].vertex[j].x * cos(B) * cos(C) - stl->facet_start[i].vertex[j].y * sin(C) * cos(B) + stl->facet_start[i].vertex[j].z * sin(B);
        ny = stl->facet_start[i].vertex[j].x * sin(C) + stl->facet_start[i].vertex[j].y * cos(C);
        nz =-stl->facet_start[i].vertex[j].x * sin(B) * cos(C) + stl->facet_start[i].vertex[j].y * sin(B) * sin(C) + stl->facet_start[i].vertex[j].z * cos(B);
        stl->facet_start[i].vertex[j].x=nx;
        stl->facet_start[i].vertex[j].y=ny;
        stl->facet_start[i].vertex[j].z=nz;
    }
  }
  
  // carefull... we are using the part size!
  stl_get_size(stl);
}

//void
//my_stl_rotate_xy(stl_file *stl, float x, float y) {
//  int i;
//  int j;
//  double nx, ny, nz;
//
//  if (stl->error) return;
//
//    for(i = 0; i < stl->stats.number_of_facets; i++) {
//    for(j = 0; j < 3; j++) {
//      my_rotate(&stl->facet_start[i].vertex[j].y,
//                 &stl->facet_start[i].vertex[j].z, x);
//      my_rotate(&stl->facet_start[i].vertex[j].z,
//                 &stl->facet_start[i].vertex[j].x, y);
//    }
//  }
//  
//  // carefull.... we are using size of the part
//  stl_get_size(stl);
//  
//}

int place3DPrinterHead(stl_file* head1, stl_file* head2, v3* location, stl_file* head1_adjusted, stl_file* head2_adjusted) {
    
    if(!head1 || !head2 || !head1_adjusted || !head2_adjusted){
        perror("Null pointer in place3DPrinterHead");
        return(-1);
    }
    
    // initialize local stl structure
    stl_initialize(head1_adjusted);
    stl_initialize(head2_adjusted);

    // make a copy of stats
    head1_adjusted->stats = head1->stats;
    head2_adjusted->stats = head2->stats;

    // allocate memory for facets
    head1_adjusted->facet_start = (stl_facet*) malloc(head1_adjusted->stats.number_of_facets * sizeof (stl_facet));
    head2_adjusted->facet_start = (stl_facet*) malloc(head2_adjusted->stats.number_of_facets * sizeof (stl_facet));

    if (head1_adjusted->facet_start == NULL || head2_adjusted->facet_start == NULL) {
        perror("place3DPrinterHead");
        return -1;
    }

    // copy facets
    memcpy(head1_adjusted->facet_start, head1->facet_start, head1->stats.number_of_facets * sizeof (stl_facet));
    memcpy(head2_adjusted->facet_start, head2->facet_start, head2->stats.number_of_facets * sizeof (stl_facet));

    stl_translate_relative(head1_adjusted, location->x, location->y, location->z + printer_char.extrusion_height/2.0f);
    stl_translate_relative(head2_adjusted, location->x, 0, location->z + printer_char.extrusion_height/2.0f);

    return 0;
}

int place3DPrinterTable(stl_file* table, stl_file* table_adjusted, float B, float C) {
    
    if(!table || !table_adjusted){
        perror("Null pointer in place3DPrinterTable");
        return -1;
    }
    // initialize local stl structure
    stl_initialize(table_adjusted);

    // make a copy of stats
    table_adjusted->stats = table->stats;

    // allocate memory for facets
    table_adjusted->facet_start = (stl_facet*) malloc(table_adjusted->stats.number_of_facets * sizeof (stl_facet));

    if (table_adjusted->facet_start == NULL) {
        perror("place3DPrinterTable");
        return -1;
    }

    // copy facets
    memcpy(table_adjusted->facet_start, table->facet_start, table->stats.number_of_facets * sizeof (stl_facet));

    my_stl_rotate(table_adjusted, B, C);

    //put the top of the table at z=0
    //stl_translate_relative(table_adjusted, 0, -(sin(DEG_TO_RAD(B)) * table->stats.max.z), -(cos(DEG_TO_RAD(C)) * table->stats.max.z));
    return 0;
}

int place3DPart(stl_file* part, stl_file* adjusted_part, float B, float C, bool newPart) {
    
    if(!part || !adjusted_part){
        perror("Null pointer in place3DPart");
        return -1;
    }
    
    // initialize local stl structure
    stl_initialize(adjusted_part);

    // make a copy of stats
    adjusted_part->stats = part->stats;

    // allocate memory for facets
    adjusted_part->facet_start = (stl_facet*) malloc(adjusted_part->stats.number_of_facets * sizeof (stl_facet));

    if (adjusted_part->facet_start == NULL) {
        perror("place3DPart");
        return -1;
    }

    // copy facets
    memcpy(adjusted_part->facet_start, part->facet_start, part->stats.number_of_facets * sizeof (stl_facet));

    // translate part so a collision is not detected if parts are side-by-side
    if(newPart)
        stl_translate_relative(adjusted_part,0,0,printer_char.extrusion_height/2.0f);
    // rotate part (we have de part offset -translation to z- already done)
    my_stl_rotate(adjusted_part, B, C);
    //stl_rotate_y(adjusted_part, B);

    return 0;
}

void computeExtremePoints(vector<v3> *H, stl_file *part){
    
    v3 lowerx=v3(part->facet_start[0].vertex[0].x,part->facet_start[0].vertex[0].y,part->facet_start[0].vertex[0].z);
    v3 upperx=lowerx, lowery=lowerx, uppery=lowerx, lowerz=lowerx, upperz=lowerx;
    
    for (int i = 0; i < part->stats.number_of_facets; i++) {
        for (int j = 0; j < 3; j++) {
            v3 point = v3(part->facet_start[i].vertex[j].x, part->facet_start[i].vertex[j].y, part->facet_start[i].vertex[j].z);
            //priority to x
            if(point.x<=part->stats.min.x){ // found with a lower x value
                if(point.x==lowerx.x){ // is it equal to what we already have?
                    if(point.y<lowerx.y || point.z<lowerx.z){ // is it better on y or z?
                        lowerx=point; // yes, accept it
                    }
                } else { // it's better, so accept it
                    lowerx=point;
                }
            }
            if(point.x>=part->stats.max.x){ // found with an upper x value
                if(point.x==upperx.x){ // is it equal to what we already have?
                    if(point.y>upperx.y || point.z>upperx.z){ // is it better on y or z?
                        upperx=point; // yes, accept it
                    }
                } else { // it's better, so accept it
                    upperx=point;
                }
            }
            
            //priority to y
            if(point.y<=part->stats.min.y){ // found with a lower y value
                if(point.y==lowery.y){ // is it equal to what we already have?
                    if(point.x<lowery.x || point.z<lowery.z){ // is it better on y or z?
                        lowery=point; // yes, accept it
                    }
                } else { // it's better, so accept it
                    lowery=point;
                }
            }
            if(point.y>=part->stats.max.y){ // found with an upper y value            
                if(point.y==uppery.y){ // is it equal to what we already have?
                    if(point.x>uppery.x || point.z>uppery.z){ // is it better on x or z?
                        uppery=point; // yes, accept it
                    }
                } else { // it's better, so accept it
                    uppery=point;
                }
            }
            
            //priority to z
            if(point.z<=part->stats.min.z){ // found with a lower z value
                if(point.z==lowerz.z){ // is it equal to what we already have?
                    if(point.y<lowerz.y || point.x<lowerz.x){ // is it better on y or x?
                        lowerz=point; // yes, accept it
                    }
                } else { // it's better, so accept it
                    lowerz=point;
                }
            }
            if(point.z>=part->stats.max.z){ // found with an upper z value            
                if(point.z==upperz.z){ // is it equal to what we already have?
                    if(point.y>upperz.y || point.x>upperz.x){ // is it better on y or z?
                        upperz=point; // yes, accept it
                    }
                } else { // it's better, so accept it
                    upperz=point;
                }
            }
        }
    }
    
    H->push_back(lowerx);
    H->push_back(upperx);
    H->push_back(lowery);
    H->push_back(uppery);
    H->push_back(lowerz);
    H->push_back(upperz);
}

void printSTL_Scad(char *filename, stl_file* file) {
    FILE *fp;

    //return;
    
    fp = fopen(filename, "w");
    if (fp){
        fprintf(fp,"module line(start, end, thickness = 1) {color(\"blue\") hull() {translate(start) sphere(thickness);translate(end) sphere(thickness);}}\n\n");
        for (int i = 0; i < file->stats.number_of_facets; i++) {
            for (int j = 0; j < 3; j++) {
                fprintf(fp, "line([%f,%f,%f],[%f,%f,%f],1);\n", file->facet_start[i].vertex[j].x,
                        file->facet_start[i].vertex[j].y, file->facet_start[i].vertex[j].z,
                        file->facet_start[i].vertex[(j + 1) % 3].x, file->facet_start[i].vertex[(j + 1) % 3].y,
                        file->facet_start[i].vertex[(j + 1) % 3].z);
            }
        }
        fclose(fp);
    }
}


void printSTL_Scad_no(char *filename, stl_file* file) {
    FILE *fp;

    return;
    
    fp = fopen(filename, "w");
    if (fp){
        fprintf(fp,"module line_no(start, end, thickness = 1) {color(\"red\") hull() {translate(start) sphere(thickness);translate(end) sphere(thickness);}}\n\n");
        for (int i = 0; i < file->stats.number_of_facets; i++) {
            for (int j = 0; j < 3; j++) {
                fprintf(fp, "line_no([%f,%f,%f],[%f,%f,%f],1);\n", file->facet_start[i].vertex[j].x,
                        file->facet_start[i].vertex[j].y, file->facet_start[i].vertex[j].z,
                        file->facet_start[i].vertex[(j + 1) % 3].x, file->facet_start[i].vertex[(j + 1) % 3].y,
                        file->facet_start[i].vertex[(j + 1) % 3].z);
            }
        }
        fclose(fp);
    }
}

int checkCollision(stl_file *table, stl_file *head1, stl_file *head2, vector<Mesh> *Meshes_part,
        vector<vector<pair<int,int>>> Lt, vector<pair<int,int>> Lc, int newPart, int newRot) {

    float B=Meshes_part->at(newPart).optimal_rotations[newRot].x;
    float C=Meshes_part->at(newPart).optimal_rotations[newRot].y;
    
    vector<vector<pair<int,int>>> Ltt = Lt;
    
    if(Lc.size())
        Ltt.push_back(Lc);
    
    //std::cout << "B=" << B << " C=" << C << endl;
            
    stl_file *tableAdjusted = new stl_file();
    if(!tableAdjusted)
        return -10;
    if(place3DPrinterTable(table, tableAdjusted, B, C)){
        stl_close(tableAdjusted);
        return -10;
    }
    
    //rotate the new part considering its rotating angles
    stl_file *newPart_adjusted = new stl_file();
    if(!newPart_adjusted)
        return -10;
    if(place3DPart(&(Meshes_part->at(newPart).stl), newPart_adjusted, B, C, true)){
        stl_close(newPart_adjusted);
        return -10;
    }

    
    vector<v3> H;
    computeExtremePoints(&H, newPart_adjusted);

    stl_file * head1_adjusted = new stl_file();
    stl_file * head2_adjusted = new stl_file();
    stl_file * adjusted_part = new stl_file();
    if(!head1_adjusted || !head2_adjusted || !adjusted_part){
        return -1;
    }

    int result = -1;
    for (auto h : H) {
        if (result != -1)break; //force foreach termination if result != 0
    
        stl_close(head1_adjusted);
        stl_close(head2_adjusted);
        if(place3DPrinterHead(head1, head2, &h, head1_adjusted, head2_adjusted)){
            result=-10;
            break;
        }

        //printf("\n//head1\n");
        //printSTL_Scad("head1_adj.scad",head1_adjusted);
        //printf("\n//head2\n");
        //printSTL_Scad("head2_adj.scad",head2_adjusted);

        for (vector<vector<pair<int,int>>>::iterator i = (Ltt).begin(); i != (Ltt).end(); i++) {
            for (vector<pair<int,int>>::iterator j = (*i).begin(); j != (*i).end(); j++) {
                stl_close(adjusted_part);
                if(place3DPart(&(Meshes_part->at((*j).first).stl), adjusted_part, B, C, false)){
                    result=-10;
                    break;
                }
                
                //printSTL_Scad("part_adj.scad", adjusted_part);
                
                if (checkSTLCollision(head1_adjusted, adjusted_part)) {
                    result = (*j).first;
                    break;
                } else if (checkSTLCollision(head2_adjusted, adjusted_part)) {
                    result = (*j).first;
                    break;
                }
            }
        }
    }


    stl_close(head1_adjusted);
    free(head1_adjusted);
    stl_close(head2_adjusted);
    free(head2_adjusted);

    stl_close(adjusted_part);
    free(adjusted_part);

    stl_close(newPart_adjusted);
    free(newPart_adjusted);
    
    stl_close(tableAdjusted);
    free(tableAdjusted);
    
    std::cout << "Collision result: " << result << endl;
    return result;
}


// check collision with head and table when building part
int checkStructureCollision(stl_file *table, stl_file *head1, stl_file *head2, stl_file *stl, v3 *Rot) {
    float B=Rot->x;
    float C=Rot->y;
    
    
    //std::cout << "B=" << B << " C=" << C << endl;
            
    stl_file *tableAdjusted = new stl_file();
    if(!tableAdjusted)
        return -1;
    if(place3DPrinterTable(table, tableAdjusted, B, C)){
        stl_close(tableAdjusted);
        return -1;
    }
    
    //printf("\n\n//adjusted table\n");
    //printSTL_Scad("adj_table.scad",tableAdjusted);
    //printSTL_Scad_no("adj_table_no.scad",&table);
    
    //rotate the new part considering its rotating angles
    stl_file *newPart_adjusted = new stl_file();
    if(!newPart_adjusted)
        return -1;
    if(place3DPart(stl, newPart_adjusted, B, C, true)){
        stl_close(newPart_adjusted);
        return -1;
    }
    //printf("\n//newPart\n");
    //printSTL_Scad_no("newPart_no.scad",newPart);
    //printSTL_Scad("newPart.scad",newPart_adjusted);

    // get all unique vertices
    vector<v3> H;
    computeExtremePoints(&H, newPart_adjusted);

    stl_file * head1_adjusted = new stl_file();
    stl_file * head2_adjusted = new stl_file();
    if(!head1_adjusted || !head2_adjusted){
        return -1;
    }
    
    int result=0;
    for (auto h : H) {
        if (result != 0)break; //force foreach termination if result != 0

        stl_close(head1_adjusted);
        stl_close(head2_adjusted);
        if(place3DPrinterHead(head1, head2, &h, head1_adjusted, head2_adjusted)){
            result=-1;
            break;
        }
        //printf("\n//head1\n");
        //printSTL_Scad("head1_adj.scad",head1_adjusted);
        //printf("\n//head2\n");
        //printSTL_Scad("head2_adj.scad",head2_adjusted);

        if(checkSTLCollision(head1_adjusted, tableAdjusted)) {
            result = 1;
            break;
        } else if (checkSTLCollision(head2_adjusted, tableAdjusted)) {
            result = 2;
            break;
        }        
    }

    stl_close(head1_adjusted);
    free(head1_adjusted);
    stl_close(head2_adjusted);
    free(head2_adjusted);
    
    stl_close(newPart_adjusted);
    free(newPart_adjusted);
    
    stl_close(tableAdjusted);
    free(tableAdjusted);
    
    //std::cout << "Structural Collision result: " << result << endl;
    return result;
}







vector<int> partsminus(vector<int> * P_bar, int part) {
    vector<int> p_hat;

    for (vector<int>::iterator it = (*P_bar).begin(); it != (*P_bar).end(); it++) {
        if ((*it) != part) {
            p_hat.push_back(*it);
        }
    }

    return p_hat;
}

vector<int> connected(set<pair<int, int>> *connections, set<int> *parts_floor,
        vector<int> *P_bar, vector<vector<pair<int,int>>> *Lt) {
    vector<int> con;

    //are there part to be built?
    if (P_bar->size() == 0) {
        //no, so return;
        return con;
    }

    if (Lt->size() == 0) {
        //no part already built, so return parts that are connected to the printer table
        for (vector<int>::iterator i = P_bar->begin(); i != P_bar->end(); i++) {
            if (find(parts_floor->begin(), parts_floor->end(), (*i)) != parts_floor->end()) {
                con.push_back(*i);
            }
        }

        return con;
    }

    
    for(vector<int>::iterator i = P_bar->begin(); i!= P_bar->end(); i++){
        // allow to build until we find a part connected to be built first
        bool isconnected=false; // not report it, until find that there is a connection with a built one
        for(set<pair<int,int>>::iterator j=connections->begin(); !isconnected && j!=connections->end();j++){
            int orig=(j->first), dest=(j->second);
            if((*i)==orig || (*i)==dest){
                // there is a connection
                if((*i)==dest){
                    orig=dest;
                    dest=(j->first);
                }
                // check that dest is on Lt
                bool found=false;
                for(vector<vector<pair<int,int>>>::iterator jj = Lt->begin(); !found && jj != Lt->end(); jj++) {
                    for (vector<pair<int,int>>::iterator ii = jj->begin(); ii != jj->end(); ii++) {
                        if((ii->first)==dest){
                            found=true;
                            break;
                        }
                    }
                }
                if(found){
                    isconnected=true;
                    break;
                }                
            }
        }
        if(isconnected)
            con.push_back(*i);
    }

    return con;
}

void print_L(vector<pair<int, int>> &Lc, vector<Mesh> *Meshes_part) {

    vector<pair<int, int>>::iterator j;
    for (j = Lc.begin(); j != Lc.end(); j++){
        std::cout << (*j).first << ", (" << RAD_TO_DEG(Meshes_part->at((*j).first).optimal_rotations[(*j).second].x) << ",";
        std::cout << RAD_TO_DEG(Meshes_part->at((*j).first).optimal_rotations[(*j).second].y) << "); ";
    }
}


void print_List(vector<vector<pair<int, int>>> &Lt, vector<Mesh> *Meshes_part) {

    vector<vector<pair<int, int>>>::iterator i;
    for (i = Lt.begin(); i != Lt.end(); i++) {
        std::cout << "\\{";
        print_L(*i,Meshes_part);
        std::cout << "\\}, ";
    }
}


void print_List_List(vector<Mesh> *Meshes_part) {
    int nBlock = 0, x;

    vector<vector<vector<pair<int, int>>>>::iterator k;
    for (k = L_f.begin(), x = 0; k != L_f.end(); k++, x++) {
        std::cout << "$\\mathcal{L}_{" << ++nBlock << "}=\\{";

        print_List(*k,Meshes_part);

        std::cout << "\\}$" << endl;
    }
}


int level(stl_file *table, stl_file *head1, stl_file *head2, 
        vector<Mesh> *Meshes_part, set<int> *parts_floor, set<pair<int, int>> *connections, 
        vector<vector<pair<int, int>>> Lt, vector<pair<int, int>> Lc, vector<int> P_bar) {

    //std::cout << endl<< endl<<"Entering Level with Lt=" << Lt.size() << " and Lc=" << Lc.size() << " and P=" << P_bar.size() << endl;

    std::cout << "Lt:";
    print_List(Lt, Meshes_part);
    
    std::cout << " Lc:";
    print_L(Lc, Meshes_part);
    std::cout << endl;
//    
//    std:cout << endl;
//    
//    if(Lc.size()>2){
//        std::cout << "Lc maior que dois" << endl;
//    }
    
    if (P_bar.size() == 0) {

        //std::cout << "Vazio" << endl;
        //we are finished
        if (Lc.size() == 0) {
            //print if we have concluded
            L_f.push_back(Lt);
            // we can find as many as we like
            if(L_f.size()<printer_char.max_Lf)
                return -1;
            else
                return -10;
        }
        
        return -1;
    }

    //find which parts in P_bar are connected with parts in Lt
    //this is previously performed and defined in P
    vector<int> P = connected(connections, parts_floor, &P_bar, &Lt);

    //no parts connected?
    if (P.size() == 0) {
        //no, so finish
        return -1;
    }

    // check that no part in P can lead to a collision, otherwise abort
    for (vector<int>::iterator it = P.begin(); it != P.end(); it++) {
        vector<vector<pair<int, int>>> Ltt=Lt;
        vector<pair<int, int>> Lcc=Lc;
        
        //P_hat is P_bar with part removed
        vector<int> P_hat = partsminus(&P_bar, *it);
                
        int nCollisions=0, partColl;
        vector<int> partColls;
        partColls.clear();
        for (int itRot =0 ; itRot< Meshes_part->at(*it).optimal_rotations.size(); itRot++) {
            std::cout << "Trying part " << (*it) << " with rotation " << itRot << " (" << Meshes_part->at(*it).optimal_rotations[itRot].x << "," << Meshes_part->at(*it).optimal_rotations[itRot].y << ")" << endl;
            if ((partColl=checkCollision(table, head1, head2, Meshes_part, Lt, Lc, (*it), itRot))<0) {
                //add part and rotation to current level of construction
                //proceed recursively

                Lcc = Lc;
                Lcc.push_back(make_pair(*it, itRot));
                
                Ltt = Lt;
                Ltt.push_back(Lcc);
                
                
                // Go recursivelly with the par in the current time slot
                if (printer_char.allow_simultaneous) {
                    //add current part to current level
                    //proceed recursively
                    int res=level(table, head1, head2, Meshes_part, parts_floor, connections, Lt, Lcc, P_hat);
                    if(res==-10) // found, so stop
                        return -10;
                    if(res>0 && *it!=res) // return until part res is found
                        return res;
                }
                
                // Go recursivelly with a new empty time slot
                vector<pair<int, int>> temp; // empty set
                temp.clear();
                 //close current level and add part to next level
                int res=level(table, head1, head2, Meshes_part, parts_floor, connections, Ltt, temp, P_hat);
                //std::cout << "Level returned" << endl;
                if(res==-10) // to stop cycle
                    return -10;
                if(res>0 && *it!=res) // return until part res is found
                    return res;
                
                break; // to need to persue with this part, if a latter collision has occoured
            } else {
                nCollisions++;
                partColls.push_back(partColl);
            }
        }
        if(nCollisions>=Meshes_part->at(*it).optimal_rotations.size()){
            // part cannot be placed due to collision
            // this part can never be included in such a Lt sequence, so abort to the part where a collision occurred
            // check what value should we return
            // first for Lc
            for (vector<pair<int,int>>::reverse_iterator j = Lc.rbegin(); j != Lc.rend(); j++) {
                for (vector<int>::iterator partIt = partColls.begin(); partIt != partColls.end(); partIt++) {
                    if (*partIt == (*j).first) {
                        std::cout << "Returning to " << (*j).first << endl;
                        return (*j).first;
                    }
                }
            }
            // now for Lt
            for (vector<vector<pair<int,int>>>::reverse_iterator i = Lt.rbegin(); i != Lt.rend(); i++) {
                for (vector<pair<int,int>>::reverse_iterator j = (*i).rbegin(); j != (*i).rend(); j++) {
                    for(vector<int>::iterator partIt=partColls.begin();partIt!=partColls.end();partIt++){
                        if(*partIt==(*j).first){
                            std::cout << "Returning to " << (*j).first << endl;
                            return (*j).first;
                        }
                    }
                }
            }
            return -1; // return one level
        }
    }
    
    // we only reach this after searching all the alternatives
    return -1;
}


