//
//  AppDelegate.h
//  useful
//
//  Created by Phil Ahrenkiel on 3/7/20.
//  Copyright © 2020 Phil Ahrenkiel. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <CoreData/CoreData.h>

@interface AppDelegate : NSObject <NSApplicationDelegate>

@property (readonly, strong) NSPersistentContainer *persistentContainer;


@end

