MODULE newComm

   use node_mod, only:  &
          nodei0,         & ! INTENT(IN)
          nodej0,         & ! INTENT(IN)
          nodeia, &
          nodeiz, &
          nodeja, &
          nodejz, &
          nodeibcon, &
          nmachs, &
          myNum

   use mem_grid, only:  &
                  ngrids
 
   implicit none

   TYPE bd
      integer :: iUpLeft
      integer :: jUpLeft
      integer :: iDownRight
      integer :: jDownRight
   end type bd
   type(bd), allocatable, dimension(:,:) :: border
   LOGICAL,allocatable,dimension(:,:) :: ThereIsSharedBorder
   
   INTEGER, parameter :: west=1,east=2,north=3,south=4

   logical,parameter :: dump=.true.
   
   type mmsg
      type(bd) :: neighbourCorners
      type(bd) :: myCorners
      integer :: msgSize
      integer :: neighbour
      TYPE(mmsg), pointer :: next => null()
   end type mmsg
   type(mmsg),pointer :: messageRoot(:),messageCurrent!
                                                      !Type: send,receive
                                                      
   TYPE msd
      type(bd) :: myCorners
      integer :: msgSize
      integer :: neighbour
      type(msd), pointer :: next => null()
   end type msd
   type(msd),pointer :: messToSendRoot(:),messToSend
   INTEGER :: message(6)

   INTEGER, allocatable, dimension(:) :: totalOfMessages,totalOfSends
   LOGICAL, allocatable, dimension(:):: ThereIsCommWith

CONTAINS

   SUBROUTINE findAndFillGhostZone(ghostzone)
      
      INTEGER,INTENT(IN) :: ghostZone
      
      INTEGER :: ng
      INTEGER :: iTULC(nMachs) !Temporary Upper left corner
      INTEGER :: jTULC(nMachs)
      INTEGER :: iTLRC(nMachs) !Temporary Lower right corner
      INTEGER :: jTLRC(nMachs)
     
!--(DMK-CCATT-INI)------------------------------------------------------------------
!      ThereIsCommWith=.false.
!      ThereIsSharedBorder=.false.
!      totalOfMessages=0
!      totalOfSends=0
!--(DMK-CCATT-FIM)------------------------------------------------------------------

      !Allocating beginning of list
      allocate(messageRoot(ngrids)) 
      allocate(messToSendRoot(ngrids)) 
      allocate(border(ngrids,4))
      allocate(ThereIsSharedBorder(ngrids,4))
      allocate(totalOfMessages(ngrids),totalOfSends(ngrids))
      allocate(ThereIsCommWith(nMachs))

!--(DMK-CCATT-INI)------------------------------------------------------------------
      ThereIsCommWith=.false.
      ThereIsSharedBorder=.false.
      totalOfMessages=0
      totalOfSends=0
!--(DMK-CCATT-FIM)------------------------------------------------------------------

      DO ng=1,ngrids
         messageCurrent=>messageRoot(ng) !Point to beggining of list
         messToSend=>messToSendRoot(ng)
         
         !Looking for borders that is shared with other processors
         !The subroutine do not identified which processor
         !This routine will fill Border
         CALL findBorders(ng,myNum,ghostZone)
         
         !Looking for corners in shared area with others threads/Processors
         !Now the processors will be identified from Border(s)
         CALL findCorners(ng,myNum,iTULC,jTULC,iTLRC,jTLRC)
         
         !Fill the list of messages to receive and send messages to
         !other processor in order to let it know about shared points
         !This routine will fill messageCurrent
         CALL fillAndSendMessages(ng,myNum,iTULC,jTULC,iTLRC,jTLRC)

         !Let's receive my neib from other processor for me know what to send
         !This routine will fill messToSend
         CALL receiveMessages(ng,myNum)
         
      END DO
      nullify(messageCurrent%next)
      nullify(messToSend%next)
      
   END SUBROUTINE findAndFillGhostZone
   
   SUBROUTINE findBorders(ng,myNum,ghostZone)
   
         INTEGER, INTENT(IN) :: ng,myNum,ghostZone
         INTEGER :: gz(ngrids,4)
         !Verifica a existência de vizinhos
         !nodeibcon(proc,grid)=0000
         !                     ||||
         !                     |||+-> west
         !                     ||+--> East
         !                     |+---> North
         !                     +----> South
         !Ajusta os cantos superiores e inferiores de cada borda para a transferencia
         !Borda oeste
         IF(iand(nodeibcon(myNum,ng),1)==0) THEN
            gz(ng,west)=ghostZone 
            border(ng,west)%iUpLeft=nodeia(myNum,ng)+nodei0(myNum,ng)-gz(ng,west)
            border(ng,west)%jUpLeft=nodeja(myNum,ng)+nodej0(myNum,ng)
            border(ng,west)%iDownRight=nodeia(myNum,ng)+nodei0(myNum,ng)-1
            border(ng,west)%jDownRight=nodejz(myNum,ng)+nodej0(myNum,ng)
            ThereIsSharedBorder(ng,west)=.true.
            IF(dump) WRITE(*,fmt='("W:",5(I3.3,1X))') myNum,border(ng,west)%iUpLeft, &
                                    border(ng,west)%jUpLeft,border(ng,west)%iDownRight, &
                                    border(ng,west)%jDownRight;CALL flush(6)
         END IF      
         !Borda leste
         IF(iand(nodeibcon(myNum,ng),2)==0) THEN
            gz(ng,east)=ghostZone
            border(ng,east)%iUpLeft=nodeiz(myNum,ng)+nodei0(myNum,ng)+1
            border(ng,east)%jUpLeft=nodeja(myNum,ng)+nodej0(myNum,ng)
            border(ng,east)%iDownRight=nodeiz(myNum,ng)+nodei0(myNum,ng)+gz(ng,east)
            border(ng,east)%jDownRight=nodejz(myNum,ng)+nodej0(myNum,ng)
            ThereIsSharedBorder(ng,east)=.true.
            IF(dump) WRITE(*,fmt='("E:",5(I3.3,1X))') myNum,border(ng,east)%iUpLeft, &
                                             border(ng,east)%jUpLeft,border(ng,east)%iDownRight, &
                                             border(ng,east)%jDownRight;CALL flush(6)
         END IF  
         !Borda norte
         IF(iand(nodeibcon(myNum,ng),4)==0) THEN
            gz(ng,north)=ghostZone 
            border(ng,north)%iUpLeft=nodeia(myNum,ng)+nodei0(myNum,ng)
            border(ng,north)%jUpLeft=nodeja(myNum,ng)+nodej0(myNum,ng)-gz(ng,north)
            border(ng,north)%iDownRight=nodeiz(myNum,ng)+nodei0(myNum,ng)
            border(ng,north)%jDownRight=nodeja(myNum,ng)+nodej0(myNum,ng)-1
            ThereIsSharedBorder(ng,north)=.true. 
            IF(dump) WRITE(*,fmt='("N:",5(I3.3,1X))') myNum,border(ng,north)%iUpLeft, &
                                             border(ng,north)%jUpLeft,border(ng,north)%iDownRight, &
                                             border(ng,north)%jDownRight;CALL flush(6)
         END IF
         !Borda sul
         IF(iand(nodeibcon(myNum,ng),8)==0) THEN
            gz(ng,south)=ghostZone
            border(ng,south)%iUpLeft=nodeia(myNum,ng)+nodei0(myNum,ng)
            border(ng,south)%jUpLeft=nodejz(myNum,ng)+nodej0(myNum,ng)+1
            border(ng,south)%iDownRight=nodeiz(myNum,ng)+nodei0(myNum,ng)
            border(ng,south)%jDownRight=nodejz(myNum,ng)+nodej0(myNum,ng)+gz(ng,south)
            ThereIsSharedBorder(ng,south)=.true.
            IF(dump) WRITE(*,fmt='("S:",5(I3.3,1X))') myNum,border(ng,south)%iUpLeft, &
                                             border(ng,south)%jUpLeft,border(ng,south)%iDownRight, &
                                             border(ng,south)%jDownRight;CALL flush(6)
         END IF

   
   end subroutine findBorders
   
   SUBROUTINE findCorners(ng,myNum,iTULC,jTULC,iTLRC,jTLRC)
   
      INTEGER, INTENT(IN) :: ng,myNum 
      INTEGER, INTENT(OUT), DIMENSION(nMachs) :: iTULC,jTULC,iTLRC,jTLRC
      
      INTEGER :: dir,iMyPos,jMyPos,neib,iNbPos,jNbPos

      iTULC=999999
      jTULC=999999
      iTLRC=0
      jTLRC=0

      Do dir=west,south !(through west,east,north,south)
         IF(.not. ThereIsSharedBorder(ng,dir)) CYCLE !There are not to do in this direction - is BC
         !Walk in all points of my ghostzone at "dir" head
         DO iMyPos=border(ng,dir)%iUpLeft,border(ng,dir)%iDownRight
            DO jMyPos=border(ng,dir)%jUpLeft,border(ng,dir)%jDownRight
               !Loop in other threads
               DO neib=1,nMachs !All processor/Threads
                  IF(neib==myNum) CYCLE !Except myself
                  !Searching my i,j point in i,j points at other threads
                  DO iNbPos=nodeia(neib,ng)+nodei0(neib,ng),nodeiz(neib,ng)+nodei0(neib,ng)
                     DO jNbPos=nodeja(neib,ng)+nodej0(neib,ng),nodejz(neib,ng)+nodej0(neib,ng)
                        IF(iMyPos==iNbPos .and. jMyPos==jNbPos) THEN !Found a shared point
                           !Store the neighbour Upper Left Corner (iTULC,jTULC)
                           !in temporary for processor neib
                           iTULC(neib)=min(iTULC(neib),iNbPos)
                           jTULC(neib)=min(jTULC(neib),jNbPos)
                           !Store the neighbour lower right Corner (iTRLC,jTRLC)
                           !in temporary for processor neib
                           iTLRC(neib)=max(iTLRC(neib),iNbPos)
                           jTLRC(neib)=max(jTLRC(neib),jNbPos)
                           !
                           !Notice that there is communication with neib processor
                           ThereIsCommWith(neib)=.true.
                        END IF
                     END DO                  
                  END DO
               END DO
            END DO
         END DO
      END DO

   END subroutine findCorners
   
   subroutine fillAndSendMessages(ng,myNum,iTULC,jTULC,iTLRC,jTLRC)
   
      include "mpif.h"
      INTEGER, INTENT(IN) :: ng,myNum 
      INTEGER, INTENT(IN), DIMENSION(nMachs) :: iTULC,jTULC,iTLRC,jTLRC
      
      INTEGER :: neib
      INTEGER :: message(6)
      INTEGER :: ierr
      INTEGER :: status(MPI_STATUS_SIZE)
      
      !Put all messages in a 
      DO neib=1,nMachs
         
         IF(ThereIsCommWith(neib)) THEN !There is comm with neib
            
            IF(dump) WRITE(*,fmt='("Comm:",6(I3.3,1X))') myNum,neib,iTULC(neib), &
                                                      jTULC(neib),iTLRC(neib),jTLRC(neib);CALL flush(6)
            !Fill message points from neighbour
            messageCurrent%neighbour=neib
            messageCurrent%neighbourCorners%iUpLeft=iTULC(neib)   -nodei0(neib,ng)
            messageCurrent%neighbourCorners%jUpLeft=jTULC(neib)   -nodej0(neib,ng)
            messageCurrent%neighbourCorners%iDownRight=iTLRC(neib)-nodei0(neib,ng)
            messageCurrent%neighbourCorners%jDownRight=jTLRC(neib)-nodej0(neib,ng)
            !messageCurrent%msgSize
            !Fill message point in my domain
            messageCurrent%myCorners%iUpLeft=iTULC(neib)   -nodei0(myNum,ng)
            messageCurrent%myCorners%jUpLeft=jTULC(neib)   -nodej0(myNum,ng)
            messageCurrent%myCorners%iDownRight=iTLRC(neib)-nodei0(myNum,ng)
            messageCurrent%myCorners%jDownRight=jTLRC(neib)-nodej0(myNum,ng)
            
            !Allocating the new node
            allocate(messageCurrent%next)
            totalOfMessages(ng)=totalOfMessages(ng)+1
            !Fill the message array to send to neib with corners
            !let it now that have comm, message(6) is 1.
            message(1)=neib
            message(2)=messageCurrent%neighbourCorners%iUpLeft
            message(3)=messageCurrent%neighbourCorners%jUpLeft
            message(4)=messageCurrent%neighbourCorners%iDownRight
            message(5)=messageCurrent%neighbourCorners%jDownRight
            message(6)=1
   
         ELSE !There is not comm with neib
            
            IF(myNum/=neib) THEN !Exclude comm with myself
               !IF not coomm with this neib processor then
               !let it now that message(6) is zero, i.e., no comm
               message(1)=neib
               message(2)=0
               message(3)=0
               message(4)=0
               message(5)=0
               message(6)=0
            END IF
   
         END IF
         
         !OK send the message to my neighbours 
         !it will charge info about the corners to be send
         IF(myNum/=neib) THEN
            CALL MPI_send(message,6, MPI_INTEGER,neib-1,10000, &
                           MPI_COMM_WORLD, ierr)
         END IF
      END DO  

   end subroutine fillAndSendMessages
   
   subroutine receiveMessages(ng,myNum)
      
      include "mpif.h"
      INTEGER, INTENT(IN) :: ng,myNum 
      
      INTEGER ::message(6)
      INTEGER :: neib
      INTEGER :: ierr
      INTEGER :: status(MPI_STATUS_SIZE)

      DO neib=1,nMachs
         IF(myNum==neib) CYCLE ! No comm with myself
         CALL MPI_recv(message,6, MPI_INTEGER,neib-1,10000, & 
               MPI_COMM_WORLD, status, ierr)
         IF(message(6)==0) CYCLE !There is no comm with this neib
         IF(dump) WRITE(*,fmt='("Recv:",6(I3.3,1X))') myNum,neib,message(1), &
                                                      message(2),message(3),message(4);CALL flush(6)
         !Fill the pointer list with my data to send
         messToSend%neighbour=neib
         messToSend%myCorners%iUpLeft=message(1)
         messToSend%myCorners%jUpLeft=message(2)
         messToSend%myCorners%iDownRight=message(3)
         messToSend%myCorners%jDownRight=message(4)
         allocate(messToSend%next)
         totalOfSends(ng)=totalOfSends(ng)+1
      END DO

   end subroutine receiveMessages

END MODULE newComm
