import React from 'react';
//import logo from './logo.svg';
import './App.css';


function ToDoForm({addTodo,incCount}){

  const [v,setValue]= React.useState("");
  let [count,setC]=React.useState(0);
  //const c=React.useState("");
  console.log("incct:", count);
  const onSubFunc= e=>{
    e.preventDefault();
    if (!v) return;
    addTodo(v);
    incCount(count++);
    setC(count);
    setValue("");
  };
  const l= 'hello';
  const formVal=(
  <form onSubmit={onSubFunc}>
  <input
    type="text"
    className="input"
    value={v}
    onChange={e=>setValue(e.target.value)}

  />
  <div>{v}</div>
  <div>{l}</div>
</form>
);

  return(
    formVal
  )

}


function App() {
  const [todos, setTodos] = React.useState([
    { data: " " ,
    isDone:false},
  ]);
  const addTodo = text => {
    const newTodos = [...todos, { text }];
    setTodos(newTodos);
  };
  const completeTodo = index => {
    const newTodos = [...todos];
    newTodos[index].isCompleted = true;
    setTodos(newTodos);
  };

  const removeTodo = index => {
    let newTodos = [...todos];
    //let idx=newTodos.indexOf(newTodos[index]);
    //if (idx !== -1) {
    //  newTodos.splice(idx, 1);
    //}
    newTodos=newTodos.filter(element=>element!==newTodos[index])

    setTodos(newTodos);
  };
  let [clickCount, setCount] = React.useState(0);
  const incCount =i=> {
    clickCount=i+1;
    setCount(i+1); 
    console.log(i,clickCount)
  };
  if (clickCount!==0)
  {
  return (
    <div className="app">
      <div className="todo-list">
        {todos.map((todo, index) => (
          <Todo
            key={index}
            index={index}
            todo={todo}
            completeTodo={completeTodo}
            removeTodo={removeTodo}
          />
        ))}
        <ToDoForm addTodo={addTodo} incCount={incCount} />
      </div>
      <div>{clickCount}</div>
    </div>
  );
  }
  else{
    setTodos([]);
    return(
    <div className="app">
    <div className="todo-list">
      <ToDoForm addTodo={addTodo} incCount={incCount} />
    </div>
    <div>{clickCount}</div>
  </div>
    );
  }
}

function Todo({ todo,index,completeTodo,removeTodo }) {
  return (
    <div className="todo"
    style={{ textDecoration: todo.isCompleted ? "line-through" : "" }}
    >
      {todo.text}
      <div className="buttons">
        <button onClick={() => completeTodo(index)}>Done</button>
      </div>
      <div className="buttons2" onClick={() => removeTodo(index)}>
        <img src="delete-img.png"  width="20" height="20" alt="delete text" className="deleteimg"/>
      </div>
    </div>
  );
};


export default App;
