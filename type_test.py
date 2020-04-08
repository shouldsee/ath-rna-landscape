from datetime import datetime
from typing import List,Dict
from pydantic import BaseModel

class User(BaseModel):
    id: int
    name = 'John Doe'
    signup_ts: datetime = None
    friends: List[int] = []
    @classmethod
    def from_web(cls,data):
    	data['friends'] = data['friend_list'].split()
    	return cls.parse_obj(data)

external_data = {
    'id': '123',
    'other':'extra',
    'signup_ts': '2019-06-01 12:22',
    'friends': [1, 2, '3']
}


web_data = {
    'id': '123',
    'other':'extra',
    'signup_ts': '2019-06-01 12:22',
    'friend_list': '1 2 3',
}


user = User(**external_data)
user = User.from_web(web_data)
print(user)