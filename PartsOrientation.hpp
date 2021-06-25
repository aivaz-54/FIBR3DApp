/** 
 * @file PartsOrientation.hpp
 * @author Bruna
 * @date Created on 22 de Maio de 2019, 15:29
 * @brief Methods to calculate the order to print each part of the object.
 * @details During the procedure there are methods to check if there are collisions between the part and the printer in order to create just the sequences that the printer is able to physically print.
 */

#ifndef PARTSORIENTATION_HPP
#define PARTSORIENTATION_HPP

#include "Inspection.hpp"
#include "admesh/stl.h"
#include "FIBR3Dapp.hpp"

using namespace std;

/**
 * @brief Place the two parts that compose the printer head
 * @param[in] head1 a part of the printer head
 * @param[in] head2 another part of the printer head
 * @param[in] location location where the printer should put the filament
 * @param[out] head1_adjusted part adjusted to the location point
 * @param[out] head2_adjusted another part adjusted to the location point
 * @return integer with the status of the adjustment
 * @retval 0 parts were adjusted correctly
 * @retval -1 error getting stl_file information
 * @note do not forget to free head1_adjusted and head2_adjusted
 * @warning it is assumed the the point [0,0,0] of stl_file head1 is where the filament will be deposited (before adjustment)
 */
int place3DPrinterHead(stl_file* head1, stl_file* head2, v3* location, stl_file* head1_adjusted, stl_file* head2_adjusted);

/**
 * @brief Place the table
 * @param[in] table stl_file with the table
 * @param[out] table_adjusted adjusted table considering B and C parameters
 * @param[in] B rotation angle in y axis
 * @param[in] C rotation angle in z axis
 * @return integer with the status of the adjustment
 * @retval 0 parts were adjusted correctly
 * @retval -1 error getting stl_file information
 * @note do not forget to free table_adjusted
 * @warning it is assumed the the point [0,0,0] of stl_file table is the center 
 */
int place3DPrinterTable(stl_file* table, stl_file* table_adjusted, float B, float C);

/**
 * @brief rotate the part according to the defined rotation angles
 * @details this method assumes that all previous parts are also rotated
 * @details the entire piece will be printed in its initial position
 * @param[in] part stl_file with the part to be rotated
 * @param[out] adjusted_part the part of the object with the rotations applied
 * @param[in] B rotation angle along y axis
 * @param[in] C rotation angle along z axis
 * @param[in] newPart if part is built or to be built
 * @return integer with the status of the adjustment
 * @retval 0 part was adjusted correctly
 * @retval -1 error getting stl_file information
 * @note do not forget to free adjusted_part
 */
int place3DPart(stl_file* part, stl_file* adjusted_part, float B, float C, bool newPart);

/**
 * @brief Check for possible collisions between the printer, all previous printed parts and the new part to be printed
 * @param[in] table_adjusted stl_file with the table in the position to print the new part
 * @param[in] head1 part of the printer head (without any adjustment)
 * @param[in] head2 another part of the printer head (without any adjustment)
 * @param[in] newPart part of the object to be printed
 * @param[in] partsOfObject all previous printed parts without any rotation
 * @param[in] B rotation angle along y axis
 * @param[in] C rotation angle along z axis
 * @return integer with the status of the adjustment
 * @retval 0 part was adjusted correctly
 * @retval 1 collision between head1 and table
 * @retval 2 collision between head2 and table
 * @retval 3 collision between existing part and head1
 * @retval 4 collision between existing part and head2
 * @retval 5 collision between existing part and table
 */
int collide3DPrinter(stl_file* table_adjusted, stl_file* head1, stl_file* head2, stl_file* newPart,
        vector<stl_file*> partsOfObject, float B, float C);

/**
 * @brief Print the stl_file in SCad format
 * @param[in] file object to be printed
 */
void printSTL_Scad(char *filename, stl_file* file);
void printSTL_Scad_no(char *filename, stl_file* file);

/**
 * @brief Creates a new set from the original without a specific element
 * @param[in] P_bar list of parts not assigned to any level
 * @param[in] part id of the part to be removed
 * @return a new set without the selected id of the part
 */
vector<int> partsminus(vector<int> * P_bar, int part);

/**
 * @brief Get a set with all parts that may be connected to the current building sequence
 * @param[in] connections set with all possible connections between all parts that compose the entire stl_file
 * @param[in] parts_floor set of all parts that are placed on the floor, i.e., directly in the table
 * @param[in] P_bar list of parts not assigned to any level
 * @param[in] Lt temporary list of building sequences
 * @return set of integers with the id (position) of the part in the initial vector
 */
vector<int> connected(vector<pair<int, int>> *connections, vector<int> * parts_floor,
        vector<int> * P_bar, vector<vector<int>> *Lt);


/**
 * @brief Calculate all possible printing sequences considering all parts, connections between parts and possible rotations for each part.
 * @param[in] table stl_file of the printing table
 * @param[in] head1 stl_file of the printing head part responsible for depositing the material
 * @param[in] head2 stl_file of the printing head part that supports head1
 * @param[in,out] L_f set of all possible printing sequences
 * @param[in] Meshes_part vector containing the Mesh of all parts that made the object
 * @param[in] parts_floor set identifying the index in Meshes_part of each part that are on the floor
 * @param[in] connections set of pairs defining the printing hierarchy of the parts
 * @param[in,out] Lt temporary list of building sequences
 * @param[in,out] Lc current list (the current printing level parts)
 * @param[in,out] P_bar list of parts not assigned to any level
 * @retval return -1 to persue, -10 to stop, and positive to return to part 
 */
int level(stl_file* table, stl_file* head1, stl_file* head2,
        vector<Mesh> * Meshes_part, set<int> * parts_floor, set<pair<int, int>> *connections,
        vector<vector<pair<int,int>>> Lt, vector<pair<int,int>> Lc, vector<int> P_bar);


/**
 * @brief Check for possible collisions between table, already deposited parts and the new part to be inserted.
 * @details This method ensures that there is no collision between the table and the head. Then checks for collision between parts that had already be inserted. Finaly, tests if positioning the new part will not cause any collision.
 * @param[in] table stl_file of the printing table
 * @param[in] head1 stl_file of the printing head part responsible for depositing the material
 * @param[in] head2 stl_file of the printing head part that supports head1
 * @param[in] Meshes_part vector containing the Mesh of all parts that made the object
 * @param[in] Lc current list (the current printing level parts)
 * @param[in] newPart the position in the vector Meshes_part of the part to be inserted
 * @param[in] Lt temporary list of building sequences
 * @return integer with the status of the collision
 * @retval negative for no collition
 * @retval non negative for the part number that lead to collision
 */
int checkCollision(stl_file *table, stl_file *head1, stl_file *head2, vector<Mesh> *Meshes_part,
        vector<vector<pair<int,int>>> Lt, vector<pair<int,int>> Lc, int newPart, int newRot);

/* @retval zero for no collision */
int checkStructureCollision(stl_file *table, stl_file *head1, stl_file *head2,
        stl_file *stl, v3 *Rot);

void computeExtremePoints(vector<v3> *H, stl_file *part);

void print_List_List(vector<Mesh> *Meshes_part);
#endif /* PARTSORIENTATION_HPP */

